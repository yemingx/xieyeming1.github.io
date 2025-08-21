python=/mnt/software/anaconda3/envs/ML/bin/python

date=250819
email_raw="545350945@qq.com"
email=${email_raw//@/_} 
email=${email//./_} 

mkdir -p tmp_${date}
$python gmail_fetch.py -i ${email_raw} -o tmp_${date}/${email}_${date}

time $python doubao_extract_path.py -i tmp_${date}/${email}_${date} -o tmp_${date}/${email}_${date}_fetched.path
time $python doubao_extract_csv.py -i tmp_${date}/${email}_${date} -o tmp_${date}/${email}_${date}_fetched.csv   

# md to excel
# /mnt/software/anaconda3/envs/ML/bin/python md_to_excel.py -i tmp_250819/wangjn_bgi.com_250710.md -o tmp_250819/wangjn_bgi.com_250710.xls
