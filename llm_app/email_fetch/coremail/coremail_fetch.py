import imaplib
import email
from email.header import decode_header
import argparse
from operator import itemgetter
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input', required=True, help = '')
parser.add_argument('-o', '--output', required=True, help = '')
args = parser.parse_args()

# 连接配置
IMAP_SERVER = "mailc.genomics.cn"  # 替换为实际服务器地址
USERNAME = "XXXXXXXXXX"
PASSWORD = "XXXXXXXXXX"  # 邮箱密码或授权码
TARGET_SENDER = args.input  # 目标发件人 wangjn@bgi.com

try:
    # 建立SSL加密连接
    mail = imaplib.IMAP4_SSL(IMAP_SERVER, 993)
    mail.login(USERNAME, PASSWORD)
    mail.select("INBOX")  # 选择收件箱[1,7](@ref)

    # 搜索来自目标发件人的邮件（按时间倒序取最新）
    status, msg_ids = mail.search(None, f'FROM "{TARGET_SENDER}"')
    if status != "OK" or not msg_ids[0]:
        print("未找到符合条件的邮件")
        mail.logout()
        exit()
    
    # 获取邮件ID列表（最新邮件在末尾）# 取最后一封（最新）[7](@ref)
    email_ids = msg_ids[0].split()
    latest_email_id = email_ids[-1]  

    # 获取邮件原始内容
    status, msg_data = mail.fetch(latest_email_id, "(RFC822)")
    raw_email = msg_data[0][1]
    email_message = email.message_from_bytes(raw_email)

except Exception as e:
    print(f"操作失败: {e}")
    exit()

# 解析主题（处理中文编码）
print(decode_header(email_message["Subject"]))
subject, encoding = decode_header(email_message["Subject"])[0]
if isinstance(subject, bytes):
    subject = subject.decode(encoding if encoding else "utf-8", "ignore")

# 获取发件人、日期
sender = email_message.get("From")
date = email_message.get("Date")

# 解析正文（优先纯文本部分）
body = ""
if email_message.is_multipart():
    for part in email_message.walk():
        if part.get_content_type() == "text/plain":
            body = part.get_payload(decode=True).decode("utf-8", "ignore")
            break
else:
    body = email_message.get_payload(decode=True).decode("utf-8", "ignore")

# 输出结果
print(f"【最新邮件】来自 {TARGET_SENDER}")
print(f"主题: {subject}")
print(f"日期: {date}")
print(f"正文:\n{body[:]}...")  # 截取前200字符 if needed

# write TARGET_SENDER subject date body to args.output
with open(args.output, 'w') as o:
    tmp = [TARGET_SENDER, subject, date, body]
    if not tmp[0]:
        print("未找到目标发件人")
        mail.logout()
        exit()
    o.write("\t".join(tmp) + "\n")


# 关闭连接
mail.logout()

