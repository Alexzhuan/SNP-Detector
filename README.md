# SNP Detector

## Requirements

环境要求：python3

第三方包要求：

- matplotlib
- numpy
- tqdm
- pysam
- pandas

可以通过下面命令安装所需要的包

```shell
pip install -r requirements.txt

# 如下载速度忙，可临时换源
pip install -r requirements.txt -i https://mirrors.aliyun.com/pypi/simple/ 
```

## Data Preparation

待分析的数据需要为排序好的二进制文件 `<file_name>.bam`，以及相应的索引文件`<file_name>.bam.bai` （需要与BAM文件在同一个文件夹下）。SAM文件可以通过下面命令转换为二进制BAM文件：

```shell
samtools view -S -b <file_name>.sam > <file_name>.bam
```

BAM文件排序可以通过下面命令：

```shell
 samtools sort -m 1G -o <file_name>.sorted.bam -@ 3 <file_name>.bam
```

BAM文件的索引生成可以通过以下命令：

```shell
samtools index <file_name>.sorted.bam <file_name>.sorted.bam.bai
```

另外还需要SNP文件 `<file_name>.txt` ，其字段要求如下：

| [column 1] | [column 2] | [column 3] | [column 4] | [column 5] |
| :--------: | :--------: | :--------: | :--------: | :--------: |
|     1      |    chrI    |    2914    |     1      |    G/A     |
|     2      |    chrI    |    2970    |     1      |    T/C     |
|     3      |    chrI    |    3328    |     1      |    T/C     |

其中 `[column 2]` 为染色体， `[column 3]` 为SNP位点，`[column5]` 为参考基因组的碱基对。

## Running Commands

准备好数据后，就可以运行脚本了，运行命令如下：

```shell
usage: python recom_cal.py --bam_file <file_name>.bam --snp_file <file_name>.txt --mode <mode string>

optional arguments:
    --batch BATCH_NUM                   设置Batch数量，默认为8
    --process PROCESS_NUM               设置进程数量，默认为8 
    --visual                            可视化结果
    --test                              测试模式
    --lb LEFT_BOUNDARY                  测试模式下，设置处理数据的左侧索引（即从第几条开始），默认是39836
    --rb RIGHT_BOUNDARY                 测试模式下，设置处理数据的右侧索引（即从第几条结束），默认是39840
```

其中，`--bam_file` 参数为设定bam文件的路径，可以选择绝对路径 `<user_root>/<project_name>/<file_name>.bam`，可以选择相对路径  `<file_name>.bam` ，要求是脚本和数据在同一文件夹下。同样，`--snp_file` 参数为设定SNP文件的路径。

`--mode` 参数为设定发生切换的条件限制，其输入需符合格式 `<mode_string,adjacent_l_k,adjacent_r_k>`, 如 "AB,0,0" 的含义为，SNP位点A，其若来自参考序列，则其下一位SNP位点B则需要来自另外一套序列，才视为发生切换，而其中 '0,0' 分别为表示需要A位点前面 `adjacent_l_k` 位SNP位点与A保持一致，需要B位点前面 adjacent_r_k 位SNP位点与A保持一致，如"AAAABBBB,3,3"，则是A前三位需要与A保持一致，B后三位需要与B保持一致。

**Note:** 考虑到存在当前待考察的SNP位点对中，可能存在在 reads 的 query sequence 上缺失的情况，例如，部分SNP位点为 `[2010, 2020, 2030]` ，待考察SNP位点对 `(2010, 2020)` ，但是位点 `2020` 在 query sequence 上缺失，则需要考虑上一位SNP位点或下一位SNP位点，构成新的SNP位点对 `(2010, 2030)` 。

为了加速运行速度，考虑了多进程运行，将数据分配给多个子进程处理。默认情况下，子进程数为8，可以通过 `--process ` 参数来设置子进程数量。关于数据分配的方式，原先采用一条染色体的数据作为一个task，但因为染色体的SNP位点数量大小不一，考虑可能影响运行速度。后来采用整个数据均分为 n 个 batch ，一个 batch 的数据作为一个task，可通过 `--batch `参数来设置 batch 的个数。

**Note:** 如果硬件设备允许，可以调大 `--process` 和 `--batch` 的参数值，最好 `--process` 和 `--batch` 的数值相等。关于 `--process` `--batch` 的设置，相应运行速度有如下参考：

| The way of allocation | Process | Batch |  Mode  |     Time     |
| :-------------------: | :-----: | :---: | :----: | :----------: |
|     Average Batch     |    4    |   4   | AAABBB |  0h 19m 12s  |
|     Average Batch     |    8    |   8   | AAABBB |  0h 10m 7s   |
|     Average Batch     |   16    |  16   | AAABBB | **0h 7m 2s** |

如果需要对分析结果进行可视化，可以通过设置 --visual 参数来输出可视化的结果，这里是默认输出 `result_[mode]_visual.pdf`  矢量图文件。

Note: 此处不建议在测试模式下使用 --visual 参数，因为目前可视化还无法根据数据的情况来动态调整可视化图，只有SNP数据中是包含16条染色体数据的时候才能运行正常。

为了方便测试，加入了 `--test` 参数，并且设置 `--lb` 和 `--rb` 参数，根据SNP文件数据的索引，确定需要处理的数据区间 `[LEFT_BOUNDARY, RIGHT_BOUNDARY)` ，这里右侧是开区间，默认是` LEFT_BOUNDARY=39836` , `RIGHT_BOUNDARY=39840` 。

## Results

最后的结果会输出两个文件，生成的结果文件为 `result_<mode_string>_read.xlsx` 和 `result_<mode_string>_read.xlsx` ，**保存在和BAM文件同一目录下**。

其中 `result_<mode_string>_read.xlsx` 中为符合条件的SNP位点及其相应 reads 信息，数据字段如下：

|               read_id                | chrom | match_snp1 | match_snp2 | match_snp_idx1 | match_snp_idx2 | first_read_snp | last_read_snp | read_snp                                                     |                           snp_str                            |
| :----------------------------------: | :---: | :--------: | :--------: | :------------: | :------------: | :------------: | :-----------: | :----------------------------------------------------------- | :----------------------------------------------------------: |
| f83613b9-2033-4101-83cf-37041390760b | chrI  |    3328    |    3359    |       2        |       3        |      2914      |     11730     | 2914,2970,3328,3359,3462,3688,3723,3942,4034,4067,4119,4145,4182,4324,4508,4597,4801,4818,4819,4831,4836,4859,4949,4993,5047,5218,5237,5451,5463,5466,5487,5569,5576,5627,5646,5718,5741,5756,5763,5814,6001,6031,6091,6153,6201,6356,6365,6556,6962,7076,7167,7268,7403,7646,10806,11045,11107,11610,11730 | BBB[3328:3359]AAABBBB\*BB-BBAAA\*-\*BBB\*BBBB\*BBABBBBBABBAB\*BBABBBBBB\*BA\*A |
| 55eaac57-57da-4352-9f9b-6a3e60a10316 | chrI  |    3462    |    3688    |       4        |       5        |      2914      |     11730     | 2914,2970,3328,3359,3462,3688,3723,3942,4034,4067,4119,4145,4182,4324,4508,4597,4801,4818,4819,4831,4836,4859,4949,4993,5047,5218,5237,5451,5463,5466,5487,5569,5576,5627,5646,5718,5741,5756,5763,5814,6001,6031,6091,6153,6201,6356,6365,6556,6962,7076,7167,7268,7403,7646,10806,11045,11107,11610,11730 | \*AA\*A[3462:3688]BBBBBBBBBBBABBBBBABBBBBBBBBABBBBBBBBB\*BBBBBBBBBBB\*BBBB |
| 0ceb0436-179d-4a3c-8976-1591f780f020 | chrI  |    3688    |    3723    |       5        |       6        |      2914      |     7646      | 2914,2970,3328,3359,3462,3688,3723,3942,4034,4067,4119,4145,4182,4324,4508,4597,4801,4818,4819,4831,4836,4859,4949,4993,5047,5218,5237,5451,5463,5466,5487,5569,5576,5627,5646,5718,5741,5756,5763,5814,6001,6031,6091,6153,6201,6356,6365,6556,6962,7076,7167,7268,7403,7646 | \*ABAAA[3688:3723]BBBA\*BBBBBBBBBAA-BBBBBBBBBBAA\*\*BBBB\*\*BBBBBBBBB\*B |

- **read_id**： read的ID；
- **chrom**： 染色体编号；
- **match_snp1**： 考察的SNP位点对的前一个SNP位点；
- **match_snp2**： 考察的SNP位点对的后一个SNP位点；
- **match_snp_idx1** 和 **match_snp_idx2** ：考察的SNP位点对在 SNP文件 `<file_name>.txt` 的下标，即对应 `[column 1]`；
- **first_read_snp** ：reads 的 query sequence  中第一个 SNP位点；
- **last_read_snp** ：reads 的 query sequence 中最后一个 SNP位点；
- **read_snp** ： reads 的 query sequence 中包含所有 SNP 位点，其中 -1 表示该SNP位点在 query sequence 上缺失；
- **snp_str**：reads 的 query sequence 上所有SNP位点的参考序列信息，`A` 表示当前碱基来自父系参考序列，`B` 表示当前碱基母系参考序列，`-` 表示当前碱基既不属于父系也不属于母系，`*` 表示当前碱基对应的 SNP 位点在 query sequence 上缺失，并在发生切换的位置进行了标注，即 `[SNP1:SNP2]`。

`result_<mode_string>_read_count.xlsx` 中为统计信息，数据字段如下：

|     chrom_snp_pair      | count |
| :---------------------: | :---: |
| chrXIII,[716031,716317] |  19   |
|  chrXVI,[25632,26200]   |  19   |
|  chrII,[488463,488565]  |  18   |
|  chrXI,[225909,226659]  |  17   |
| chrXIII,[720435,721068] |  16   |
|  chrII,[482219,482508]  |  15   |
| chrXII,[321387,322259]  |  10   |

- **chrom_snp_pair** ： 由 `chrom`  、`match_snp1` 和 `match_snp2` 构成的标识；
- **count**：当前 SNP 位点对发生切换的reads总数。

