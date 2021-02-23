# v4

import os
import pysam
import tqdm
import time
import pandas as pd
import datetime
import argparse
from multiprocessing import Pool
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument('--bam_file', type=str, help='bam file')
parser.add_argument('--snp_file', type=str, help='snp_file')
parser.add_argument('--mode', type=str, default='AAAAABBBBB,4,4', help='single mode')
parser.add_argument('--batch', type=int, default=8, help='batch')
parser.add_argument('--process', type=int, default=8, help='process')
parser.add_argument('--visual', action='store_true', help='visualize result')
# parser.add_argument('--read_detail', action='store_true', help='match read detail')
parser.add_argument('--test', action='store_true', help='test mode')
parser.add_argument('--lb', type=int, default=39836, help='left boundary')
parser.add_argument('--rb', type=int, default=39840, help='right boundary')

args = parser.parse_args()

snp_file = pd.read_table(args.snp_file, header=None)

contigs = snp_file[1].values
snp_list = snp_file[2].values
ref_snp_bases = [(item.split('/')[0], item.split('/')[1]) for item in snp_file[4].values]


def split_chrom():

    chrom_list = snp_file[1].unique()
    chrom_group = {}
    for chrom in chrom_list:
        filtered_snp_list = list((snp_file[snp_file[1] == chrom][2]).values)
        filtered_pos_list = list((snp_file[snp_file[1] == chrom][0] - 1).values)
        chrom_group[chrom] = (filtered_snp_list, filtered_pos_list)

    return chrom_group


chrom_snp_dict = split_chrom()
mode, adjacent_l_k, adjacent_r_k = args.mode.split(',')[0], int(args.mode.split(',')[1]), int(args.mode.split(',')[2])


def binary_search(nums, data):
    length = len(nums)
    first = 0
    last = length - 1
    while first <= last:
        mid = (last + first) // 2
        if nums[mid] is None:
            next_pos = mid
            for i in range(mid - 1, -1, -1):
                if nums[i] is None:
                    continue
                else:
                    if nums[i] > data:
                        next_pos = i
                    elif nums[i] == data:
                        return i
                    break

            if next_pos != mid:
                last = next_pos
            else:
                for i in range(mid + 1, length):
                    if nums[i] is None:
                        continue
                    else:
                        if nums[i] < data:
                            next_pos = i
                        elif nums[i] == data:
                            return i
                        break
                if next_pos != mid:
                    first = next_pos
                else:
                    return -1
        else:
            if nums[mid] > data:
                last = mid - 1
            elif nums[mid] < data:
                first = mid + 1
            else:
                return mid
    return -1

def is_same_chromosome(cur_pos, next_pos):
    cur_contig = contigs[cur_pos]

    flag = True
    for i in range(1, adjacent_l_k + 1):
        if cur_contig != contigs[cur_pos - i]:
            flag = False
            break
    for i in range(0, adjacent_r_k + 1):
        if cur_contig != contigs[next_pos + i]:
            flag = False
            break

    return flag

def parse_mode(mode, adjacent_l, adjacent_r):
    mode_len = len(mode)

    mode_l_seq = [0 if mode[i] == 'A' else 1 for i in range(adjacent_l)]
    mode_r_seq = [0 if mode[(adjacent_l + 1) + i] == 'A' else 1 for i in range(adjacent_r)]

    return mode_l_seq, mode_r_seq


def combine_expected_adjacent_base(expected_base_pair, ref_adjacent_l_snp_base1_list,
                                   ref_adjacent_l_snp_base2_list,
                                   ref_adjacent_r_snp_base1_list, ref_adjacent_r_snp_base2_list):

    expected_mode_map = {}

    expected_mode_map['adjacent_l_base_map'] = {expected_base_pair[0]: [], expected_base_pair[1]: []}
    expected_mode_map['adjacent_r_base_map'] = {expected_base_pair[0]: [], expected_base_pair[1]: []}

    mode_l_seq, mode_r_seq = parse_mode(mode, adjacent_l_k, adjacent_r_k)
        
    for idx, item in enumerate(mode_l_seq):
        if not item:
            expected_mode_map['adjacent_l_base_map'][expected_base_pair[0]].append(ref_adjacent_l_snp_base1_list[idx])
            expected_mode_map['adjacent_l_base_map'][expected_base_pair[1]].append(ref_adjacent_l_snp_base2_list[idx])
        else:
            expected_mode_map['adjacent_l_base_map'][expected_base_pair[0]].append(ref_adjacent_l_snp_base2_list[idx])
            expected_mode_map['adjacent_l_base_map'][expected_base_pair[1]].append(ref_adjacent_l_snp_base1_list[idx])

    for idx, item in enumerate(mode_r_seq):
        if item:
            expected_mode_map['adjacent_r_base_map'][expected_base_pair[0]].append(ref_adjacent_r_snp_base2_list[idx])
            expected_mode_map['adjacent_r_base_map'][expected_base_pair[1]].append(ref_adjacent_r_snp_base1_list[idx])
        else:
            expected_mode_map['adjacent_r_base_map'][expected_base_pair[0]].append(ref_adjacent_r_snp_base1_list[idx])
            expected_mode_map['adjacent_r_base_map'][expected_base_pair[1]].append(ref_adjacent_r_snp_base2_list[idx])

    return expected_mode_map


class Result:

    def __init__(self):
        self.match_read = {}
        self.result_count_dict = {}
        self.result_list_dict = {}

        self.result_count_dict['recom_count'] = 0
        self.result_count_dict['recom_count_AB'] = 0
        self.result_count_dict['recom_count_BA'] = 0
        self.result_count_dict['reads'] = []

        self.result_list_dict['recom_count'] = []
        self.result_list_dict['recom_count_AB'] = []
        self.result_list_dict['recom_count_BA'] = []
        self.result_list_dict['recom_rate'] = []
        self.result_list_dict['reads'] = []

        self.match_read['read_id'] = []
        self.match_read['read_snp'] = []
        self.match_read['snp_str'] = []
        self.match_read['chrom'] = []
        self.match_read['match_snp1'] = []
        self.match_read['match_snp2'] = []
        self.match_read['match_snp_idx1'] = []
        self.match_read['match_snp_idx2'] = []
        self.match_read['first_read_snp'] = []
        self.match_read['last_read_snp'] = []

    def add(self, snp_str, read_snp, read_id, contig, cur_snp, next_snp, cur_snp_idx, next_snp_idx, first_read_snp, last_read_snp, orientation='AB'):
        self.result_count_dict['recom_count'] += 1

        if orientation == 'AB':
            self.result_count_dict['recom_count_AB'] += 1
        elif orientation == 'BA':
            self.result_count_dict['recom_count_BA'] += 1

        self.result_count_dict['reads'].append(read_id)

        self.match_read['read_id'].append(read_id)
        self.match_read['read_snp'].append(read_snp)
        self.match_read['snp_str'].append(snp_str)
        self.match_read['chrom'].append(contig)
        self.match_read['match_snp1'].append(cur_snp+1)
        self.match_read['match_snp2'].append(next_snp+1)
        self.match_read['match_snp_idx1'].append(cur_snp_idx)
        self.match_read['match_snp_idx2'].append(next_snp_idx)
        self.match_read['first_read_snp'].append(first_read_snp)
        self.match_read['last_read_snp'].append(last_read_snp)

    def update(self, read_filtered_count):
        self.result_list_dict['recom_count'].append(self.result_count_dict['recom_count'])
        self.result_list_dict['recom_count_AB'].append(self.result_count_dict['recom_count_AB'])
        self.result_list_dict['recom_count_BA'].append(self.result_count_dict['recom_count_BA'])
        self.result_list_dict['recom_rate'].append(self.result_count_dict['recom_count'] / read_filtered_count)
        self.result_list_dict['reads'].append(",".join(self.result_count_dict['reads']))

        self.clear()

    def clear(self):
        self.result_count_dict['recom_count'] = 0
        self.result_count_dict['recom_count_AB'] = 0
        self.result_count_dict['recom_count_BA'] = 0
        self.result_count_dict['reads'] = []


def get_read_boundary(cur_snp, next_snp, adjacent_l_snp_list, adjacent_r_snp_list):
    left_boundary = cur_snp
    right_boundary = next_snp

    if adjacent_l_snp_list[0] < left_boundary:
        left_boundary = adjacent_l_snp_list[0]

    if adjacent_r_snp_list[-1] > right_boundary:
        right_boundary = adjacent_r_snp_list[-1]

    return left_boundary, right_boundary


def modify_l_snp_(ref_pos, cur_snp_idx, read_start, cur_contig):
    revised_l_snp_pos = -1
    revised_l_snp_idx = cur_snp_idx - 1
    flag = 0
    while revised_l_snp_idx >= adjacent_l_k and (snp_list[revised_l_snp_idx] - 1) >= read_start and contigs[revised_l_snp_idx] == cur_contig:
        revised_l_snp_pos = binary_search(ref_pos, (snp_list[revised_l_snp_idx] - 1))
        if revised_l_snp_pos != -1:
            flag = 1
            break
        else:
            revised_l_snp_idx -= 1
    if flag == 0:
        return -1, -1, -1
    else:
        return snp_list[revised_l_snp_idx] - 1, revised_l_snp_pos, revised_l_snp_idx
    
def modify_r_snp_(ref_pos, cur_snp_idx, read_end, cur_contig):
    revised_r_snp_pos = -1
    revised_r_snp_idx = cur_snp_idx + 1
    flag = 0
    while revised_r_snp_idx < len(snp_list) - adjacent_r_k and (snp_list[revised_r_snp_idx] - 1) <= read_end and contigs[revised_r_snp_idx] == cur_contig:
        revised_r_snp_pos = binary_search(ref_pos, (snp_list[revised_r_snp_idx] - 1))
        if revised_r_snp_pos != -1:
            flag = 1
            break
        else:
            revised_r_snp_idx += 1
    if flag == 0:
        return -1, -1, -1
    else:
        return snp_list[revised_r_snp_idx] - 1, revised_r_snp_pos, revised_r_snp_idx
    
def modify_snp_pair(ref_pos, cur_snp, next_snp, cur_snp_idx, read_start, read_end):
    cur_snp_query_pos = binary_search(ref_pos, cur_snp)
    next_snp_query_pos = binary_search(ref_pos, next_snp)
    next_snp_idx = cur_snp_idx + 1
    cur_contig = contigs[cur_snp_idx]
    
    if cur_snp_query_pos == -1:
        cur_snp, cur_snp_query_pos, cur_snp_idx = modify_l_snp_(ref_pos, cur_snp_idx, read_start, cur_contig)
    if next_snp_query_pos == -1:
        next_snp, next_snp_query_pos, next_snp_idx = modify_r_snp_(ref_pos, next_snp_idx, read_end, cur_contig)

    if cur_snp_query_pos == -1 or next_snp_query_pos == -1:
        return None
    
    ref_cur_snp_base1, ref_cur_snp_base2 = ref_snp_bases[cur_snp_idx]
    ref_next_snp_base1, ref_next_snp_base2 = ref_snp_bases[next_snp_idx]
    
    adjacent_l_snp_list = []
    adjacent_r_snp_list = []
    ref_adjacent_l_snp_base1_list = []
    ref_adjacent_l_snp_base2_list = []
    ref_adjacent_r_snp_base1_list = []
    ref_adjacent_r_snp_base2_list = []
    adjacent_l_snp_query_pos = []
    adjacent_r_snp_query_pos = []
    adjacent_l_snp_idx_list = []
    adjacent_r_snp_idx_list = []
        
    tmp_adjacent_l_snp_idx = cur_snp_idx - 1
    tmp_adjacent_r_snp_idx = next_snp_idx + 1
    
    for i in range(adjacent_l_k):
        if tmp_adjacent_l_snp_idx < 0:
            return None
        else:
            if (snp_list[tmp_adjacent_l_snp_idx] - 1) < read_start or contigs[tmp_adjacent_l_snp_idx] != cur_contig:
                return None
        
        adjacent_l_cur_snp = snp_list[tmp_adjacent_l_snp_idx] - 1
        adjacent_l_cur_snp_query_pos = binary_search(ref_pos, adjacent_l_cur_snp)
        
        if adjacent_l_cur_snp_query_pos == -1:
            adjacent_l_cur_snp, adjacent_l_cur_snp_query_pos, tmp_adjacent_l_snp_idx = modify_l_snp_(ref_pos, tmp_adjacent_l_snp_idx, read_start, cur_contig)
        if adjacent_l_cur_snp_query_pos == -1:
            return None
        else:
            adjacent_l_snp_list.append(adjacent_l_cur_snp)
            ref_adjacent_l_snp_base1_list.append(ref_snp_bases[tmp_adjacent_l_snp_idx][0])
            ref_adjacent_l_snp_base2_list.append(ref_snp_bases[tmp_adjacent_l_snp_idx][1])
            adjacent_l_snp_query_pos.append(adjacent_l_cur_snp_query_pos)
            adjacent_l_snp_idx_list.append(tmp_adjacent_l_snp_idx)
            tmp_adjacent_l_snp_idx -= 1
    
    for i in range(adjacent_r_k):
        if tmp_adjacent_r_snp_idx >= len(snp_list):
            return None
        else:
            if (snp_list[tmp_adjacent_r_snp_idx] - 1) > read_end or contigs[tmp_adjacent_r_snp_idx] != cur_contig:
                return None
        
        adjacent_r_cur_snp = snp_list[tmp_adjacent_r_snp_idx] - 1
        adjacent_r_cur_snp_query_pos = binary_search(ref_pos, adjacent_r_cur_snp)
        
        if adjacent_r_cur_snp_query_pos == -1:
            adjacent_r_cur_snp, adjacent_r_cur_snp_query_pos, tmp_adjacent_r_snp_idx = modify_r_snp_(ref_pos, tmp_adjacent_r_snp_idx, read_end, cur_contig)
        if adjacent_r_cur_snp_query_pos == -1:
            return None
        else:
            adjacent_r_snp_list.append(adjacent_r_cur_snp)
            ref_adjacent_r_snp_base1_list.append(ref_snp_bases[tmp_adjacent_r_snp_idx][0])
            ref_adjacent_r_snp_base2_list.append(ref_snp_bases[tmp_adjacent_r_snp_idx][1])
            adjacent_r_snp_query_pos.append(adjacent_r_cur_snp_query_pos)
            adjacent_r_snp_idx_list.append(tmp_adjacent_r_snp_idx)
            tmp_adjacent_r_snp_idx += 1
    
    expected_base_pair = [(ref_cur_snp_base1, ref_next_snp_base2), (ref_cur_snp_base2, ref_next_snp_base1)]
    expected_mode_map = combine_expected_adjacent_base(expected_base_pair, ref_adjacent_l_snp_base1_list,
                                                        ref_adjacent_l_snp_base2_list,
                                                        ref_adjacent_r_snp_base1_list,
                                                        ref_adjacent_r_snp_base2_list)
    
    return (cur_snp, next_snp, cur_snp_idx, next_snp_idx, cur_snp_query_pos, next_snp_query_pos, adjacent_l_snp_list, adjacent_r_snp_list,
           adjacent_l_snp_query_pos, adjacent_r_snp_query_pos, adjacent_l_snp_idx_list, adjacent_r_snp_idx_list, expected_base_pair, expected_mode_map)


def is_adjacent_snp_match_ref(read, adjacent_l_snp_query_pos, adjacent_r_snp_query_pos, adjacent_l_snp_idx_list, adjacent_r_snp_idx_list):
    adjacent_l_query_bases = []
    adjacent_r_query_bases = []

    for idx, query_pos in enumerate(adjacent_l_snp_query_pos):
        if read.query_sequence[query_pos] in ref_snp_bases[adjacent_l_snp_idx_list[idx]]:
            adjacent_l_query_bases.append(read.query_sequence[query_pos])
        else:
            return False
    for idx, query_pos in enumerate(adjacent_r_snp_query_pos):
        if read.query_sequence[query_pos] in ref_snp_bases[adjacent_r_snp_idx_list[idx]]:
            adjacent_r_query_bases.append(read.query_sequence[query_pos])
        else:
            return False

    return adjacent_l_query_bases, adjacent_r_query_bases


def search_read_snp(snp_list, ref_start, ref_end):
    length = len(snp_list)
    first = 0
    last = length - 1

    snp_l_pos = None
    snp_r_pos = None
    while first < last:
        mid = (last + first) // 2
        if snp_list[mid] < ref_start:
            first = mid + 1
        else:
            last = mid
    snp_l_pos = last

    first = 0
    last = length - 1
    while first < last - 1:
        mid = (last + first) // 2
        if snp_list[mid] <= ref_end:
            first = mid
        else:
            last = mid
    snp_r_pos = first

    return snp_l_pos, snp_r_pos


def cal_adjacent_snp_match_expected(query_base_pair, expected_mode_map, adjacent_l_query_bases,
                                    adjacent_r_query_bases, result, expected_base_pair, read, ref_pos, contig, cur_snp, next_snp, cur_snp_query_pos, next_snp_query_pos, cur_snp_idx, next_snp_idx):

    if expected_mode_map['adjacent_l_base_map'][query_base_pair] == adjacent_l_query_bases and expected_mode_map['adjacent_r_base_map'][query_base_pair] == adjacent_r_query_bases:
        if query_base_pair == expected_base_pair[0]:
            orient = 'AB'
        else:
            orient = 'BA'

        ref_start = read.reference_start
        ref_end = read.reference_end

        cur_chrom_snp_list, cur_chrom_pos_list = chrom_snp_dict[contig]
        snp_l_pos, snp_r_pos = search_read_snp(cur_chrom_snp_list, ref_start, ref_end)
        read_snp_list = [snp - 1 for snp in cur_chrom_snp_list[snp_l_pos:snp_r_pos + 1]]
        read_snp_query_pos = [binary_search(ref_pos, snp) for snp in read_snp_list]

        first_read_snp = read_snp_list[0]+1
        last_read_snp = read_snp_list[-1]+1

        reads_mode_str = ""
        for idx, query_pos in enumerate(read_snp_query_pos):
            
            if query_pos == -1:
                reads_mode_str += '*'
            elif read.query_sequence[query_pos] == ref_snp_bases[cur_chrom_pos_list[snp_l_pos] + idx][0]:
                reads_mode_str += 'A'
            elif read.query_sequence[query_pos] == ref_snp_bases[cur_chrom_pos_list[snp_l_pos] + idx][1]:
                reads_mode_str += 'B'
            else:
                reads_mode_str += '-'

            if query_pos == cur_snp_query_pos:
                reads_mode_str += f'[{cur_snp+1}:{next_snp+1}]'

        read_snp = ','.join([str(snp+1) for snp in read_snp_list])

        result.add(reads_mode_str, read_snp, read.query_name, contig, cur_snp, next_snp, cur_snp_idx, next_snp_idx, first_read_snp, last_read_snp, orient)


def is_switch(batch_id, bam_file, lb=39254, rb=39894):
    print("{0}\tAnalyze Batch {1} [{2}, {3}]...".format(datetime.datetime.now(), batch_id, lb, rb-1))
    result = Result()

    base_list = []
    chrom_list = []
    read_count_list = []  # read count list

    start_time = time.time()

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for i in tqdm.tqdm(range(len(contigs))[lb: rb]):
            contig = contigs[i]

            # switch case: - ? A  A  A  B  B  B ? - or more

            # judge whether current snp and adjacent snp locate in same chromosome
            if not is_same_chromosome(i, i + 1):
                continue

            cur_snp = snp_list[i] - 1
            next_snp = snp_list[i + 1] - 1

            # adjacent left index orientation: 0 <- 1 <- 2 <- ... <- adjacent_l_k - 1 <- cur
            # adjacent right index orientation: next -> 0 -> 1 -> 2 -> ... -> adjacent_r - 1
            adjacent_l_snp_list_ = [snp - 1 for snp in snp_list[i - adjacent_l_k:i]]
            adjacent_r_snp_list_ = [snp - 1 for snp in snp_list[i + 2:i + 2 + adjacent_r_k]]

            read_count = 0

            for read in bam.fetch(contig=contig, start=cur_snp, stop=cur_snp + 1):
                
                read_count += 1
                left_boundary, right_boundary = get_read_boundary(cur_snp, next_snp, adjacent_l_snp_list_,
                                                                  adjacent_r_snp_list_)

                if read.reference_end < right_boundary or read.reference_start > left_boundary:
                    continue

                ref_pos = read.get_reference_positions(full_length=True)
                revised_res = modify_snp_pair(ref_pos, cur_snp, next_snp, i, read.reference_start, read.reference_end)
                if revised_res is not None and isinstance(read.query_sequence, str):
                    modified_cur_snp, modified_next_snp, cur_snp_idx, next_snp_idx, cur_snp_query_pos, next_snp_query_pos, adjacent_l_snp_list, adjacent_r_snp_list, adjacent_l_snp_query_pos, adjacent_r_snp_query_pos, adjacent_l_snp_idx_list, adjacent_r_snp_idx_list, expected_base_pair, expected_mode_map = revised_res
                    
                    query_cur_base = read.query_sequence[cur_snp_query_pos]
                    query_next_base = read.query_sequence[next_snp_query_pos]

                    if query_cur_base in ref_snp_bases[cur_snp_idx] and query_next_base in ref_snp_bases[next_snp_idx]:
                        tmp = is_adjacent_snp_match_ref(read, adjacent_l_snp_query_pos, adjacent_r_snp_query_pos, adjacent_l_snp_idx_list, adjacent_r_snp_idx_list)
                        if tmp:
                            adjacent_l_query_bases, adjacent_r_query_bases = tmp
                            if (query_cur_base, query_next_base) in expected_base_pair:
                                cal_adjacent_snp_match_expected((query_cur_base, query_next_base), expected_mode_map,
                                                                adjacent_l_query_bases, adjacent_r_query_bases, result, 
                                                                expected_base_pair, read, ref_pos, contig, modified_cur_snp, modified_next_snp, cur_snp_query_pos, next_snp_query_pos, cur_snp_idx, next_snp_idx)

            result.update(read_count)
            base_list.append(cur_snp + 1)
            read_count_list.append(read_count)
            chrom_list.append(contig)

    print("{0}\tAnalyze Batch {1} [{2}, {3}] done, time consuming {4}min".format(datetime.datetime.now(), batch_id, lb, rb-1, (time.time() - start_time)/60))

    return result, base_list, read_count_list, chrom_list


def save_result(result, read_df_dict):

    read_df = pd.DataFrame({'read_id': result.match_read['read_id'],
                            'chrom': result.match_read['chrom'],
                            'snp_str': result.match_read['snp_str'],
                            'read_snp': result.match_read['read_snp'],
                            'match_snp1': result.match_read['match_snp1'],
                            'match_snp2': result.match_read['match_snp2'],
                            'match_snp_idx1': result.match_read['match_snp_idx1'],
                            'match_snp_idx2': result.match_read['match_snp_idx2'],
                            'first_read_snp': result.match_read['first_read_snp'],
                            'last_read_snp': result.match_read['last_read_snp']})
    read_df_dict.append(read_df)

    return read_df_dict


def visual_result(result_df, save_file):

    filtered_df = result_df[['chrom', 'base_pos', 'recom_count', 'recom_count_AB', "recom_count_BA"]]
    chrom_list = filtered_df['chrom'].unique()

    plt.figure(figsize=(12, 45))
    plt.subplots_adjust(hspace=0.35)
    for idx, chrom in enumerate(chrom_list):
        plt.subplot(len(chrom_list), 1, idx + 1)
        tmp_df = filtered_df[filtered_df['chrom'] == chrom]
        plt.title(chrom, fontweight='bold')
        tmp2_df = tmp_df[tmp_df['recom_count'] > 0]
        plt.scatter(tmp2_df['base_pos'], tmp2_df['recom_count_AB'], s=19, marker='o', c='', edgecolors='b')
        plt.scatter(tmp2_df['base_pos'], tmp2_df['recom_count_BA'], s=19, marker='o', c='', edgecolors='r')
        plt.grid(axis='y', linestyle='--')

    plt.savefig(save_file, bbox_inches="tight")


def split_average(batch=8):

    max_adjacent_l_k, max_adjacent_r_k = adjacent_l_k, adjacent_r_k

    batch_dict = {}
    batch_size = len(snp_list) // batch

    for i in range(batch):
        min_index = i * batch_size
        max_index = (i + 1) * batch_size

        if i == 0:
            min_index = max_adjacent_l_k

        if i == batch - 1:
            max_index = len(snp_list) - max_adjacent_r_k - 1

        batch_dict[i] = (min_index, max_index)

    return batch_dict


if __name__ == '__main__':

    chrom_group = split_average(batch=args.batch)

    pool = Pool(processes=args.process)
    result_dict = {}
    read_df_dict = []

    start_time = time.time()

    if args.test:
        result, base_list, read_count_list, chrom_list = is_switch('test', args.bam_file, args.lb, args.rb, )
        read_df_dict = save_result(result, read_df_dict)

    else:
        for chrom in chrom_group.keys():
            lb, rb = chrom_group[chrom]
            # is_switch(chrom, args.bam_file, lb, rb,)
            result_dict[chrom] = pool.apply_async(is_switch, args=(chrom, args.bam_file, lb, rb, ))
        pool.close()
        pool.join()

        for key in result_dict.keys():
            result, base_list, read_count_list, chrom_list = result_dict[key].get()
            read_df_dict = save_result(result, read_df_dict)

    result_file_name = 'result_' + mode
    dir_ = os.path.dirname(args.bam_file)

    final_read_df = pd.concat(read_df_dict)
    final_read_df.drop_duplicates(subset=['read_id', 'chrom', 'match_snp1', 'match_snp2', 'match_snp_idx1', 'match_snp_idx2', 'first_read_snp', 'last_read_snp', 'read_snp'], inplace=True)
    final_read_df.to_excel(os.path.join(dir_, result_file_name + '_read.xlsx'), columns=['read_id', 'chrom', 'match_snp1', 'match_snp2', 'match_snp_idx1', 'match_snp_idx2', 'first_read_snp', 'last_read_snp', 'read_snp', 'snp_str'],
                            index=None)

    final_read_df['chrom_snp_pair'] = final_read_df[['chrom', 'match_snp1', 'match_snp2']].apply(lambda x: x[0]+','+'['+str(x[1])+','+str(x[2])+']', axis=1)
    df = final_read_df.chrom_snp_pair.value_counts().reset_index()
    df.columns = ['chrom_snp_pair', 'count']
    df.to_excel(os.path.join(dir_, result_file_name + '_read_count.xlsx'), index=None)

    if args.visual:
        visual_result(final_df, os.path.join(dir_, result_file_name + '_visual.pdf'))

    end_time = time.time() - start_time
    print("{0}".format(datetime.datetime.now()) + "\tFinished! Total time consuming: {:.0f}h {:.0f}m {:.0f}s".format(end_time // 3600, (end_time % 3600) // 60, end_time % 3600 % 60))
