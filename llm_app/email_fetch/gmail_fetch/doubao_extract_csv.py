import argparse
from operator import itemgetter
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input', required=True, help = '')
parser.add_argument('-o', '--output', required=True, help = '')
args = parser.parse_args()

# Configuration
TEST_PROMPT = open(args.input, "r").read()

sys_prompt=""" Analyze the email content below and perform the following tasks:
Extract ALL sample name and sample ID pairs from these sections:
    - "评估不同体积的捕获效果" table (Protein/Volume sections)
    - "IL4" comparison table (Blocker/Dispersion sections)
    Format as CSV with two columns: "Sample_ID" and "Sample_Name"
    Example: 41,IL10_m5_10ul_200pg
    all special characters (such as "- + / ( ) . , : ; ! ? ' " etc) in sample name should be replaced with _
    all sample name should contain no empty space

Output format:
[CSV data here]

Email content for reference:
[Insert full email text here] """


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


