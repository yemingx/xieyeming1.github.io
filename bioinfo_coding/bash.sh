### qsub -cwd -l vf=2G,num_proc=2 -P RNAProj -binding linear:2 -q bc_rd.q,bc.q -e log -o log [bash.sh]
batchfile=($(find ${fq_dir} -name "*_demultiplex_R1.fq" -type f))
sample_id_tmp=(${batchfile[@]##*/})
sample_id=(${sample_id_tmp[@]%_demultiplex_R1.fq})
echo ${sample_id[@]}
mkdir -p chunk
for i in ${sample_id[@]}; do
  if test ! -d "chunk/${i}";then
    mkdir -p chunk/${i}
  fi
done

################## table ##################
### check freq
awk '{A[$1]++}END{for(i in A)print i,A[i]}' 
awk '{A[$1]++}END{for(i in A)print i,A[i]}' |tr [:blank:] \\t|awk '{A[$2]++}END{for(i in A)print i,A[i]}'
cat FILE|awk '{A[$4]++}END{for(i in A)print i,A[i]}' |tr [:blank:] \\t|awk '{b[$1]=$2;sum=sum+$2} END{for (i in b) print i,b[i],(b[i]/sum)*100}'|tr [:blank:] \\t|sort -k1 -nr
### dedup
awk '!x[$0]++'

### selected column
awk '$2 > '0' && $4 == 'C''
awk '$2 == 'A' ; $4 == 'C'' 
awk '($2=="chr7" && $4!="chr7");($2!="chr7" && $4=="chr7")'
thres=20
awk -v t="$thres" '$2 < t'
awk -v a="$a" -v b="$b" '$1==a && $10 == b'  
### string partly match
awk '$3 ~ /snow/ { print }'     # contain certain words
lr |awk '$5 ~ /M/ { print }'|awk '{print $9}'|xargs rm
awk '{if($2~/GCCT$/)print $0}'
awk '{if($1~/^[ABCD]/)print $0}'
awk '$1!~/[a-z]/'

### sort column rank:
sort -k4,4 -k1,1 -k2,2n

### print number of columns, rows:
awk '{print NF}'
wc -l

### print nth line:
sed -n '1;3p' 
sed -n '1,3p'
sed -n '3p'
sed '3q;d'

### conditional column print
awk '{if ($3==$4) print $5"\t"$3;else print $5"\t"$3"\n"$5"\t"$4}'
awk '{if ($4 >= 0.3) {print $1"\t"$2"\t"$3"\t"0.3}else {print}}' hek_nome3h_45M_chr7_m_ratio_GpC.0_07.bedgraph
awk '{ if (($9 - $8) > 60) print }'

### column substring, use awk to print line, for column 5, print the 5th to 7th character for the string after _
awk '{split($5, arr, "_"); print $1, $2, $3, $4, substr(arr[2], 5, 3)}' filename

### left join
awk 'NR==FNR{A[$1]=$2;next}{print$0 FS (A[$1]?A[$1]:"missing")}' file1 file2
### output: common_column, file2_column, file1_column

### awk to read two files at the same time
awk 'FNR==NR { # actions for the first file } FNR!=NR { # actions for the second file }' file1 file2

################## line check ##################
### print every second line
awk 'NR%4==2'
awk '{if (NR%4==2) print length}' 

### print line length:
awk '{print length}'
expr length AGCAGCCGT

### grep line start with, end with
grep ^a*
grep ^@
grep -v ^@
grep "^>"
grep "/1$" (grep end with /1)

### grep A or B:
egrep -h "stringA|stringB" filename
grep 'AGGCAAAC\|CTAACTTA\|TACGGGCG\|GCTTTCGT'
grep 'wordA*'\''wordB' filename
top|grep -w -e 'bedmap' -e 'awk'

### grep line with given line length range
grep -x '.\{3,10\}'

### grep same line from two files, show the line in file2
grep -Ff file1 file2
### grep same string from two file, show the line in file2
grep -wf file1 file2

### grep if contain string
### grep if contain string, leading or trailing N characters
grep -Eo '.{0,10}ACATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCAAGCAGAAGACGGCATACGAGCTCTTCCGATGT.{0,10}' file.txt
### grep if contain string, at least appear N times
grep -E '(GATC.*){2,}'

###Count Occurrences of Character per Line
head sonic_100pg.clusters|awk -F 'chr7' 'BEGIN{print "Line", "\tCount"}{print NR "\t" NF-1}'
###Count Occurrences of Character per field
cut -d '|' -f 2 items.txt | tr -c -d "E\n" | cat -n
################## line edit ##################
### see line ending \r
cat -vE file
od -c file
### to remove \r at line ending
sed 's/\r$//'
sed 's/.$//'

### to remove empty line
awk NF file.txt

### trim line to length and to position:
cut -c 12 file
cut -c 2-13 file

### grep selected line
sed '2,4!d' somefile.txt

### add string to the end of the first line of every 4 lines
sed '1~4 s/$/string/g'
sed '1~4 s/_/ /2'
sed 's/$/\n\n\n/g'
awk '{print $1"\n""\n""\n"}'
sed '1~4 s/^/string/g'

sed -n '1~4p;2~4p'

cat ${R1} |sed 'N;N;N;s#\n#\t#g'|awk '{print $1"\n"$4"\n+\n"$8}'|sed '1~4 s/\./\_/' > ${sample}_raw_R1.fq
cat ${R2} |sed 'N;N;N;s#\n#\t#g'|awk '{print $1"\n"$4"\n+\n"$8}'|sed '1~4 s/\./\_/'|sed '1~4 s/\.2/\.1/' > ${sample}_raw_R2.fq

### add string to the end of every line split by tab
awk '{print $0 "\tEcoRI"}' file.txt

### remove last n character of each line
sed 's/.....$//'
echo 987654321 | rev | cut -c 4- | rev

### remove line ending
tr -d "\n\r"

### remove leading and trailing whitespace 
### and also squeeze sequences of tabs and spaces into a single space
awk '{$1=$1};1'

### merge two file line after line
paste -d $'\n' file1 file2

### merge two gz file
cat file1.gz file2.gz > allfiles.gz

### replace word:
sed 's/word1/word2/g'

### replace in certain column
awk 'BEGIN {FS=OFS="\t"} {if ($5 == 0) $5 = "+"} 1'

### replace the selected occurence:
sed 's/\t/\_/4'
sed 's/\S*value$/test/' # replace value of last column, \S for is for non-whitespace, $ is for the end of line.

### convert fastq to fasta:
sed -n '1~4s/^@/>/p;2~4p' 

### convert 4 line fastq to one line fastq:
sed 'N;N;N; s/\n/ /g'
sed 'N;N;N;s#\n#\t#g'
### convert one line fastq back to 4 line fastq
sed 's/\ /\n/g'
### convert multiple line to one line
tr '\n' '\t'
tr -d '\n'

### convert tab to newline
sed -e 'y/\t/\n/'

### use variable in sed, use double quote
### Variables inside ' don't get substituted in Bash
### Single quotes won't interpolate anything, but double quotes will.
a='foo'         
echo 'baz' | sed 's/baz/'"$a"'/g'

### remove first n lines
sed -e '1,3d'
tail -n +4
awk 'NR > 3 { print }'

### remove selected lines
sed -e '10d;11d;13d;14d' read_dist.tmp1

### random lines from a file
shuf -n N input > output

### remove duplicated line
awk '!x[$0]++'

### add line number starting from 1
awk '{print $3"\t"$5"\t"1+i++}'

################## math ##################
### quant pct
awk '{b[$1]=$2;sum=sum+$2} END{for (i in b) print i,b[i],(b[i]/sum)*100}'

### sum of column
awk '{s+=$1} END {print s}' 

### calculate mean:
awk '{ total += $1 } END { print total/NR }'
mean=$(printf "%s\n" "${numbers[@]}" | awk '{ sum += $1 } END { printf "%.3f", sum / NR }')

### Calculate the standard deviation
sd=$(printf "%s\n" "${numbers[@]}" | awk -v mean="$mean" '{ sum += ($1 - mean) ^ 2 } END { printf "%.3f", sqrt(sum / NR) }')

### column calc
awk '{print $1*100"\t"$2"\t"$3}'
awk '{print $3 - $2}'

### basic int calc
echo $((2**15))

### int division
printf %.3f\\n "$((1000 * 2/15))e-3"
awk "BEGIN {printf \"%.2f\n\", 100/3}"
awk '{printf "%d\n", int((log($1)/log(2)-10)/0.125) }'
cat hic_pro_metrics.xls|cut -f1,6,11|grep -v sample|awk '{printf "%s\t%.3f\n",$1,$3/$2}'

ratio=($(awk "BEGIN {printf \"%.2f\n\", $hg19_mt /4 /$E }"))

### float calc
echo 2*0.15 | bc
echo 1121/134343|bc -l
echo 'print $(( 10/3.3 ))' | zsh

### column divide to decimal
awk '{printf "%.2f%s",$7/$8,"\n"}'

### min in array
min=0
for num in "${sample_id[@]}";do if ((num < min));then min=$num;fi;done
echo $min

### max in array
max=0
for num in "${sample_id[@]}";do if ((num > max));then max=$num;fi;done
echo $max

### ceiling
echo "2" | awk '{print int($1+0.999999)}'

### floor
awk '{print int(i++/2)}'

################## string ##################
### string length
a='AGCCCTAACCCTCCG'
echo ${#a}

### remove sam tag
sed 's/\tOQ\:Z\:[^\t]*//'
sed 's/\tXN\:[^\t]*//'

### replace blank space with tab
tr [:blank:] \\t
tr -s "\t" < read_dist.tsv

### remove repeat string
tr -s ' '

### insert character to string
sed 's/./&,/4'

### Remove character before/ after a pattern
$ a='hello:world'
$ b=${a%:*}
$ echo "$b"
hello
$ a='hello:world:of:tomorrow'
$ echo "${a%:*}"
output: hello:world:of
$ echo "${a%%:*}"
output: hello
$ echo "${a#*:}"
output: world:of:tomorrow
$ echo "${a##*:}"
output: tomorrow

### match files
ls chunk_{1..10} #  brace expansion is used to generate a sequence of numbers

################## array ##################
seqfile=($(find ${_07_bt2_whitelist_dir} -name "*_read.id" -type f))
seqfile=($(ls ${_07_bowtie2_whitelist_dir}/read_id/*_read.id))
sample_id_tmp=(${seqfile[@]##*/})	# keep string after last /
sample_id=(${sample_id_tmp[@]%_read.id})	# trim suffix

${X[@]}                                        the usual whole array
${X[@]:index:length}                           slice from index to index+length-1 inclusive
${X[@]::length}                                slice from 0 to length-1 inclusive
${X[@]:index}                                  slice from index to end of array inclusive
X=( "${X[@]}" )                                compact X
X=( "${X[@]::$INDEX}" "${X[@]:$((INDEX+1))}" ) remove element at INDEX and compact array

Y=( 000 111 222 333 )
echo "${Y[@]:2}"
output: 222 333

echo "${Y[@]:1:3}"
output: 111 222 333

echo "${Y[@]::1}"
output: 000

### convert string to array
echo "abcdefg" | fold -w1
echo "abcdefg" | grep -o .

################## file ##################
### find directory or file
find . -type d -name "name*"
find . -type f -name "name*"

### change all the permission to 755
find . -type d -exec chmod 755 {} \;

### find multiple file extension
find . -type f \( -iname \*.fq -o -iname \*.fastq \)

### change permission recursively


## remove zero byte fastq
rm -f $(lr |awk '$5==0'|awk '{print $9}')
find . -name '*.fq' -size 0 -print0 | xargs -0 rm
find . -name 'file*' -size 0 -delete

## remove files in subfolders
find . -type f -name '*sh.e*'|xargs rm
find . -name "*.bam*" -type f|grep -v hg19|xargs rm

## find specific batch of sample and do bash operation
find . -type f -name *sanger_to_trinity.sam|xargs wc -l
### if src/ only contains files:
find src/ ! -name Default.png -exec cp -t dest/ {} +
### If src/ has sub-directories, this omits them, but does copy files inside of them:
find src/ -type f ! -name Default.png -exec cp -t dest/ {} +
find hichew_*/ -type f -name *vecs.tsv -exec cp -t eigen_ins_sum/ {} +
find helaBMpECTCFGly016_chr7_* -name "*.ctcf.xls.gz" -exec cp -t . {} +
### If src/ has sub-directories, this does not recurse into them:
find src/ -type f -maxdepth 1 ! -name Default.png -exec cp -t dest/ {} +

### ls files ONLY
ls -p | grep -v /

### copy whole directory files and folders
cp -R <source_folder>/* <destination_folder>

### compress file
gzip -c sample.fq > sample.fq.gz
### decompress file
gzip -dc < sample.fq.gz > /somewhere/sample.fq
zcat *fastq.gz > merged.fq
unzip experimental_queries-master.zip

### compress folder or file
tar -czvf project.tar.gz /path/to/project
tar -cjf ARCHIVE_NAME.tar.bz2 [FILE_NAMES] 

### decompress all files from .tar or .tar.gz
tar -xf project.tar.gz -C /somewhere/

### decompress all files from bzip2
tar -xjf project.tar.bz2 -C /somewhere/

### view tar file content
tar -tvf archive.tar
tar -ztvf archive.tar.gz
tar -jtvf archive.tar.bz2
tar -Jtvf archive.tar.xz

### split file by line or byte
split -l 10 input.txt /somewhere/prefix
split -b 10 input.txt /somewhere/prefix

### search "text" in current directory and all the files inside 
grep -R "COMBINED_SA.SUBGEM.bsv.gff.GEMLINE_OUT.tsv.BE2_FRAGNUM" .

################## system ##################
### remove in batch
qstat |awk '{print $1}' |xargs qdel
qstat |grep fastp_bt2| awk '{print $1}' |xargs qdel
lr HelaBM_CTCF_pC0320_gpc_*/HelaBM_CTCF_pC0320_gpc.10.100.*X.ctcf.*|awk '$7==13'|awk '{print $9}'|xargs rm

### remove files over 1M
find . -type f -size +1M -exec rm {} +

### self check file number and directory size
ls -lR /zfswh1/BC_RD_P0/P19Z12200N0089_xym|grep "^-"|wc -l
du -sh *

### check disk size and quota
df -h 
lfs quota -gh bc_tsrd /opt/zfswh1
lfs quota -gh bc_tsrd /zfssz5/
qstat -F vf,cpu,mt -q bc_rd.q

### check run memory and real time
/usr/bin/time -f "%E %M" *<command>*
/usr/bin/time -v command  #Maximum resident set size
### check gpu usage
gpustat -i

### clean screen: 
clear

### create symbolic link
ln -s my_file.txt my_link.txt
### override
ln -sf my_file.txt my_link.txt
### unlink
rm my_link.txt

### check difference between two files
diff file1 file2

### sort out top files
ls -Slh chunk/|grep fq|head -1200|awk '{print $9}'> L01_scSPRITE.top1200
files=($(cat L01_scSPRITE.top1200))
mkdir -p top1200
target_dir="top1200"
for file in "${files[@]}"; do cp "chunk/$file" "$target_dir"; done

### nohup pid check
nohup my_command > my.log 2>&1 &
echo $! > save_pid.txt
nohup sh bash.sh &

### check pid working dir
pwdx [PID]

### print stdout and stderr to terminal and to file
() 2>&1 | tee file

### How do I save (append) output to an existing output file?
echo "Today's date is" >> foo.txt
date >> foo.txt

### check cpu gpu ram usage, process
number of cpu: cat /proc/cpuinfo | grep processor | wc -l
CPU: top
ps -ef |grep nohup 
RAM: free -m
nvidia-smi -l 1

### ls time order
ls -lcth

### grep certain process
ps aux|grep 
### kill signal
kill -9 [PID]

### chmod -R 755 path
### before -rw-r--r--   after -rwxr-xr-x


### Run script1.sh in the background
./script1.sh &
### Run script2.sh in the background
./script2.sh &
### Run script3.sh in the background
./script3.sh &
### Wait for all scripts to finish
wait

### install GCC
sudo yum group install "Development Tools"
sudo yum install man-pages

### mount nas drive
mount -t nfs XXXXXXXXXXX

################## remote login ##################
ssh user@host ls -l /some/directory
ssh phoenixnap@185.52.53.222

su chenweitian
chenwt12345qwert

################## rsync ##################
### download
rsync -av xieyeming1@10.0.0.14:[source] [destination] 

### upload
rsync -av [source] xieyeming1@10.0.0.14:[destination] 

rsync [OPTION]... [USER@]HOST:SRC [DEST]

myVar=${myVar:0:1}${myVar:2:1}
echo "ff9f75d4e7bda792fca1f30fc03a5303  package.deb" | md5sum -c –

################## ftp ##################
ftp ftp.cngb.org
pass
ls
delete
put file

HOST='ftp.cngb.org'
USER='ngb_00616'
PASSWD='dGTcjHLX@0'

ftp -p -n $HOST <<END_SCRIPT
quote USER $USER
quote PASS $PASSWD
put CvipI_4T1_NOM_PC_7_5_0929.rd.fastq.gz
END_SCRIPT

### https://stackoverflow.com/questions/45523098/bash-script-to-upload-files-to-a-ftp-server
### https://www.cs.colostate.edu/helpdocs/ftp.html
############### .bashrc ###############
### set default python version:
alias python='/ifswh1/BC_PUB/biosoft/pipeline/Package/Python-3.5.1/python'
### reload bashrc settings
source ~/.bashrc
. ~/.bashrc
### set terminal showing path:
PS1='\u@\h: \w $ '
### custom commands
alias lr="ls -lhtcr"
alias lx="ls -lhtcrX"

############### write module ###############
### change your PATH variable to add /home/USER/miniconda3/bin, equal to PATH=/home/USER/miniconda3/bin:$PATH (in bash or sh)
###%Module1.0
prepend-path PATH /home/USER/miniconda3/bin
### To change your LD_LIBRARY_PATH  variable to add /home/USER/geos/lib
###%Module1.0
prepend-path LD_LIBRARY_PATH /home/USER/geos/lib
###%Module1.0
module use --append /hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/module_file/
prepend-path PATH /hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/bioinfo/longranger-2.2.2/

### load module
module load /hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/module_file/longranger/2.2.2
### check loaded module
module list
############### md5 ###############
check md5
linux: md5sum
mac: md5
find . -type f -exec md5 {} + > checklist.chk

############### vim ###############
:i [enter]    —insert—
:x [return]    save and exit
:q! [return]    exit without saving

ctrl + b: previous page
ctrl + f: next page

############### var ###############
var='xxx'
echo ${var}
arr=(a b c)
echo ${arr[@]}
for i in ${arr[@]}; do echo ${i};done
for i in ${arr[@]};do
  echo ${i}
done

############### increment ###############
i=0
until [ $i -gt 3 ]
do
  echo i: $i
  ((i=i+1)) # or ((i+=1)) or ((i++))
done 

############### boolean ###############
a=2462620
b=2462620
if [ "$a" -eq "$b" ];then
  echo "They're equal";
fi
echo $(( a != b ))
-eq  ==  equal
-ne  !=  not equal
-lt  less than
-le  less than or equal
-gt  greater than
-ge  greater than or equal
if [ -f "path/to/file/FILE" ] && [ ! -f "path/to/file/FILE" ];then
  echo "";
fi
if test ! -d "./hisat2_out"; then
  mkdir ./hisat2_out
fi
############### create guest user ###############
###check all user id
cut -d: -f1 /etc/passwd
getent passwd
###create guest user
sudo adduser guest_weizhi
sudo passwd -d guest_weizhi
###remove user
sudo userdel usr_id
###change password
sudo passwd usr_id

############### environment ###############
### set $PATH variable
export PATH=$PATH:/place/with/the/file
export PATH=~/opt/bin:$PATH
echo $PATH    

### create tmp
export TMPDIR=tmp
### clean sys files
find /tmp -ctime +10 -exec rm -rf {} +   
sudo apt autoremove && sudo apt autoclean
sudo du -xh --max-depth=1 /var
sudo du -xh --max-depth=1 /var/log
### Note, change +30 to the number of days you want to keep.
sudo find /var/log -mtime +30 -type f -delete

### qsub -cwd -V -l vf=1G,num_proc=1 -P RNAProj -q bc_rd.q,bc.q -e log -o log [bash.sh]
-cwd execute job from current working directory
-l resource request list
-V export all environmental variables from terminal to job
-q submit to a specific queue
-P project name for the job

### cuda version
nvcc --version

################## task specific ##################
for ((i=1; i<=10000; i++)); do   echo "$i" >> numbers.txt; done

while IFS= read -r line_number; do  sed -n "${line_number}p" numbers.txt;done < chr7_gene.index

### add header to file
awk 'BEGIN{printf "Sr No\tName\tSub\tMarks\n"} {print}' marks.txt

### grep not aligned reads map, convert to fasta
grep -v ^@ *.sam|awk '$2 == "4"'|cut -f1,10|sed -e 'y/\t/\n/'|sed 's/V/\>V/g'

### trim fix sequence
zcat ${read_w_barcode}|sed 'N;N;N; s/\n/ /g'|cut -c -47,73-|cut -c -161,187-|sed 's/\ /\n/g' \
> ${_07_trim_tag_L02_dir}/L02_R1.fq.trim_fix_seq

zcat ${read_no_barcode}|cut -c -150 > ${_07_trim_tag_L02_dir}/L02_R2.fq.trim_fix_seq

### for loop
for i in *; do echo "cat ${i} |awk '{A[$3]++}END{for(i in A)print i,A[i]}'";done > check_dup_align.txt
for i in $(seq 1 $END); do echo $i; done

counter=1
for item in "${items[@]}"; do
    echo "Item $counter: $item"
    ((counter++))  
done

seqfile=($(find /research/xieyeming1/proj_2022/fy_heteroChrom_20220706/depth/ -name "*depth" -type f))
sample_id_tmp=(${seqfile[@]##*/})
sample_id=(${sample_id_tmp[@]%.depth})
echo ${sample_id[@]}

for i in ${sample_id[@]}; do
  sample="${i}"
  if test ! -d "${sample}";then
  if [ ! -f ${sample}/chunk/${i}/${sample}_${i}_merged.bed_seg ];then
    echo "${i}"
  fi
done


### down sampling
all_cov=(16 32 64 128 256 512 1024 2048 4096 8192)
sample_id=($(ls ${test_report_dir}|grep -v Z))
echo ${sample_id[@]}
for cov in ${all_cov[@]}; do
  for j in ${sample_id[@]}; do
    echo "qsub -cwd -l vf=${ram}G,num_proc=${cpu} -P RNAProj -q bc_rd.q,bc.q -e ${proj_dir}/script/log -o ${proj_dir}/script/log ${proj_dir}/script/_02_cov_trinity_rsem_single.sh ${j} ${cov}"
    qsub -cwd -l vf=${ram}G,num_proc=${cpu} -P RNAProj -q bc_rd.q,bc.q -e ${proj_dir}/script/log -o ${proj_dir}/script/log ${proj_dir}/script/_02_cov_trinity_rsem_single.sh ${j} ${cov}
  done
done

### downsample report pair id
up_limit=3500000
step=9
arrVar=()
for i in $(seq $step -1 0);do
  denominator=$((2**$i))
  echo "$denominator"
  sampled_reads=($(awk "BEGIN {printf \"%.0f\n\", ${up_limit}/$denominator}"))
  arrVar+=($sampled_reads)
done
echo ${arrVar[@]}

### copy all certain files
find ../ -type f -name "*ctcf.xls.gz" -exec cp -t . {} +

### freq and pct
cat BD_ONT_400K.fq.cis|awk 'NR%4==1'|awk '{A[$1]++}END{for(i in A)print i,A[i]}'|tr [:blank:] \\t|awk '{A[$2]++}END{for(i in A)print i,A[i]}'
cat FILE|awk '{A[$4]++}END{for(i in A)print i,A[i]}' |tr [:blank:] \\t|awk '{b[$1]=$2;sum=sum+$2} END{for (i in b) print i,b[i],(b[i]/sum)*100}'|tr [:blank:] \\t|sort -k1 -nr

head -10000 dedup.bed.barcode_merged.hic_tmp1|awk '{A[$1]++}END{for(i in A)print i,A[i]}' |tr [:blank:] \\t|awk '{A[$2]++}END{for(i in A)print i,A[i]}'|tr [:blank:] \\t|\
awk '{b[$1]=$2;sum=sum+$2} END{for (i in b) print i,b[i],(b[i]/sum)*100}'|tr [:blank:] \\t|sort -k1 -nr

### write 0-leading line numbers
for (( num=1; num<=133; num++ ));do foo=$(printf "%03d" $num) && echo "${foo}";done

### fastq length
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' fastq
awk '{if(NR%2==0) {count++; bases += length} } END{print bases/count}' fasta

### awk if
head HiC_merge.csv|sed 's/\_/ /g'|grep -v chrM |awk '{if($6 - $3 > 0) print $0}'

### generate header
echo "chr segment_start segment_end read_id strand segment_num GpC_methy GpC_unmethy CpG_methy CpG_unmethy"|tr [:blank:] \\t > ${sample}_single_mol_methy_seg.xls

### calculate base number
zcat fastq.gz |awk 'NR%4==2'|awk '{print length($0)}'|awk '{s+=$1} END {print s}' 

### gtf to tss bed
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/gencode.vM19.annotation.gtf.gz
zcat gencode.vM19.annotation.gtf.gz | awk 'OFS="\t" {if ($3=="transcript") {if ($7 == "+") {print $1,$4-1,$4,$12,".",$7} else {print $1,$5-1,$5,$12,".",$7}}}' | tr -d '";' | sort -k1,1V -k2,2n > gencode.vM19.annotation.tss.bed

### seqtk pipe
${seqtk} trimfq -b 150 ${fq_dir}/V350187741_L01_read_2.fq.gz | ${seqtk} trimfq -e 2 -> V350187741_L01.barcode.fq
