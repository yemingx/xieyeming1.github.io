mkdir -p tmp
/mnt/software/anaconda3/envs/ML/bin/python coremail_fetch.py -i wangjn_bgi.com -o tmp/wangjn_bgi.com_250710

/mnt/software/anaconda3/envs/ML/bin/python deepseek_extract_path.py -i tmp/wangjn_bgi.com_250710 -o tmp/wangjn_bgi.com_250710.md

# md to excel
/mnt/software/anaconda3/envs/ML/bin/python md_to_excel.py -i tmp/wangjn_bgi.com_250710.md -o tmp/wangjn_bgi.com_250710.xls
