import requests
from requests.exceptions import JSONDecodeError
from bs4 import BeautifulSoup
import json
import re
import time
import argparse
from operator import itemgetter
parser = argparse.ArgumentParser(description='')
parser.add_argument('-m', '--max_articles', required=True, help = '')
parser.add_argument('-o', '--output', required=True, help = '')
args = parser.parse_args()
max_articles_ = int(args.max_articles)

# 大模型API设置（这里以OpenAI为例）
LLM_API_URL = "https://api.deepseek.com/chat/completions"
LLM_API_KEY = "XXXXXXXXXXXXXXX"  # 替换为你的API Key
# 模拟浏览器访问的请求头   https://zhuanlan.zhihu.com/p/714173074
headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
    'Cookie': 'XXXXXXXXXXXXXXX'  # 替换为实际登录后的Cookie（关键）
}

def get_wechat_articles(keywords, max_articles):
    """爬取微信公众号文章列表"""
    base_url = 'https://mp.weixin.qq.com'
    # search_url = f'{base_url}/cgi-bin/searchbiz?'
    url = "https://mp.weixin.qq.com/cgi-bin/appmsg"

    # 获取token接口（需要已登录Cookie）
    # token_url = f'{base_url}/cgi-bin/appmsg?'
    # token_params = {'action': 'list', 'lang': 'zh_CN'}
    # token_res = requests.get(token_url, params=token_params, headers=headers)
    # token = re.search(r'token=(\d+)', token_res.url).group(1)
    token='XXXXXXXXXXXXXXX'
    # 构造关键词搜索请求
    articles = []
    params = {
            'action': 'list_ex',
            'fakeid': 'XXXXXXXXXXXXXXX',
            "lang": "zh_CN",
            "f": "json",
            "ajax": "1",
            'query': keywords,
            'begin': '0',
            'count': '100',
            'token': token,
            "type": "9",
    }
    res = requests.get(url, params=params, headers=headers)
    # print(res.text)
    # biz_data = json.loads(res.text)
    # fakeid = biz_data['list'][0]['fakeid']
    
    # 获取文章列表
    page = 0
    while len(articles) < max_articles:
        appmsg_params = {
            'action': 'list_ex',
            'fakeid': 'XXXXXXXXXXXXXXX',
            "lang": "zh_CN",
            "f": "json",
            "ajax": "1",
            'query': keywords,
            'begin': page * 5,
            'count': '100',
            'token': token,
            "type": "9",
        }
        res = requests.get(url, params=appmsg_params, headers=headers)
        data = json.loads(res.text)
        # print(data)
        if not data.get('app_msg_list'):
            break
            
        for item in data['app_msg_list']:
            if any(kw in item['title'] or kw in item['digest'] for kw in keywords.split()):
                articles.append({
                    'title': item['title'],
                    'link': item['link'],
                    'publish_date': time.strftime("%Y-%m-%d", time.localtime(item['create_time']))
                })
        
        page += 1
        time.sleep(1)  # 防止请求过快
    # print(articles)
    return articles[:max_articles]

def extract_journal_with_llm(content):
    """使用大模型API提取发表期刊"""
    prompt = f"从以下学术报道中提取发表的论文名称和论文链接，用竖线'|'分隔（若文中无期刊信息则返回‘无’）：\n{content[:3000]}"  # 截断长文本
    # 使用requests库发送HTTP请求获取JSON响应
    url = LLM_API_URL
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {LLM_API_KEY}"
    }
    payload = {
        "model": "XXXXXXXXXXXXXXX",
        "messages": [
            {"role": "system", "content": "你是一个专业的学术信息提取助手"},
            {"role": "user", "content": prompt}
        ],
        "stream": False
    }
    
    try:
        response = requests.post(url, headers=headers, json=payload)
        response.raise_for_status()  # 检查HTTP错误状态码
        result = response.json()
    except requests.exceptions.RequestException as e:
        print(f"API请求错误: {e}")
        return ""
    except json.JSONDecodeError:
        print(f"JSON解析错误，响应内容: {response.text}")
        return ""
    
    # Validate response structure
    if not isinstance(result, dict) or 'choices' not in result or not isinstance(result['choices'], list) or len(result['choices']) == 0:
        print("Unexpected response structure: missing 'choices' array")
        # print(f"API Response content: {result}")
        return ""
    choice = result['choices'][0]
    if not isinstance(choice, dict) or 'message' not in choice or not isinstance(choice['message'], dict) or 'content' not in choice['message']:
        print("Unexpected response structure: missing 'message' content")
        return ""
    
    journal = choice['message']['content'].strip()
    # replace % with \t
    journal = journal.replace('|', '\t')
    return journal

def process_articles(articles):
    """处理文章并提取信息"""
    results = []
    for article in articles:
        # 获取文章详情
        res = requests.get(article['link'], headers=headers)
        soup = BeautifulSoup(res.text, 'html.parser')
        content = soup.find(class_='rich_media_content').text.strip() if soup.find(class_='rich_media_content') else ""
        
        # 使用大模型提取期刊
        journal = extract_journal_with_llm(content)
        
        results.append({
            "标题": article['title'],
            # "原文链接": article['link'],
            "发表日期": article['publish_date'],
            "发表期刊": journal
        })
        
        time.sleep(2)  # API调用限制
    
    return results

# 主执行
if __name__ == "__main__":
    keywords_list = ["受体", "感受器", "胁迫", "抗性", "适应"]
    for keywords in keywords_list:
        print("开始爬取文章列表...")
        articles = get_wechat_articles(keywords, max_articles=max_articles_)  # 爬取5篇文章
        
        print("\n开始提取文章信息...")
        results = process_articles(articles)
        
        print("\n结果：")
        for i, res in enumerate(results):
            # output to tsv
            with open(args.output, 'a') as f:
                f.write(f"{keywords}\t{res['标题']}\t{res['发表日期']}\t{res['发表期刊']}\n")
            # print(f"关键词：{keywords}")
            # print(f"{i+1}. [{res['发表日期']}] {res['标题']}")
            # print(f"   期刊: {res['发表期刊']}")

            # print(f"   链接: {res['原文链接']}\n")