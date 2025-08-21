import argparse
from operator import itemgetter
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input', required=True, help = '')
parser.add_argument('-o', '--output', required=True, help = '')
args = parser.parse_args()

# Configuration
TEST_PROMPT = open(args.input, "r").read()

sys_prompt=""" 你是一个邮件解析助手。请根据以下规则处理生物信息学交付邮件：
- **定位基础路径**：搜索 "数据路径：" 后的地理位置（如"武汉"），提取后续的绝对路径，作为基础路径。
- **定位项目路径**：根据用户指定 Lane 描述（如"L03的数据是瓮博项目的数据"），在邮件中搜索该 Lane（如"L03"）对应的 "PRJ" （例如 L02对应 P25Z11900N0021_ARAuunC) 和 在"L03 "后的 "WH"（例如 L02 对应 'WHB5EXONPEP00100816'），
                组合为  PRJ/WH 格式，例如 `P25Z11900N0021_ARAuunC/WHB5EXONPEP00100816`，作为项目路径。
- **生成完整路径**：将路径按以下格式组合：  
   `/基础路径/Zebra/项目路径/`  
   – 确保斜杠正确，末尾带斜杠。
- **输出要求**：仅返回最终路径，无额外解释。
 """

from openai import OpenAI
client = OpenAI(api_key="XXXXXXXXXXXXXXX", base_url="XXXXXXXXXX")
response = client.chat.completions.create(
    model="XXXXXXXXXX",
    messages=[
        {"role": "system", "content": sys_prompt},
        {"role": "user", "content": TEST_PROMPT},
    ],
    stream=False
)

print(response.choices[0].message.content)

# write TARGET_SENDER subject date body to args.output
with open(args.output, 'w') as o:
    o.write(response.choices[0].message.content)


