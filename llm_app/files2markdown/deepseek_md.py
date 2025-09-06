import os
import requests
from datetime import datetime
import pandas as pd
import argparse
from operator import itemgetter
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input', required=True, help = '')
parser.add_argument('-o', '--output', required=True, help = '')
args = parser.parse_args()

DEEPSEEK_API_URL = "https://api.deepseek.com/v1/chat/completions"
DEEPSEEK_API_KEY = "XXXXXXXX"  # Replace with your actual API key

def generate_markdown_with_deepseek():
    # Define headers at the start of the function
    headers = {
        "Authorization": f"Bearer {DEEPSEEK_API_KEY}",
        "Content-Type": "application/json"
    }
    
    # Find files in directories
    txt_files = sorted(f for f in os.listdir(args.input+'/files') if f.endswith('.txt'))
    tsv_files = sorted(f for f in os.listdir(args.input+'/files') if f.endswith('.tsv')) 
    png_files = sorted(f for f in os.listdir(args.input+'/files') if f.endswith('.png'))
    scripts_files = sorted(f for f in os.listdir(args.input+'/scripts') if f.endswith(('.py', '.R', '.sh')))

    # Prepare content
    content_parts = []
    
    # Process text files
    if txt_files:
        content_parts.append("## Text Content\n\n")
        for txt_file in txt_files:
            with open(os.path.join(args.input, 'files', txt_file), 'r') as f:
                content_parts.append(f"### {txt_file}\n```\n{f.read()}\n```\n\n")
    
    # Process TSV files as tables
    if tsv_files:
        content_parts.append("## Data Tables\n\n")
        for tsv_file in tsv_files:
            df = pd.read_csv(os.path.join(args.input, 'files', tsv_file), sep='\t')
            content_parts.append(f"### {tsv_file}\n\n{df.to_markdown(index=False)}\n\n")
    
    # Process PNG files
    if png_files:
        content_parts.append("## Visualizations\n\n")
        for i, png_file in enumerate(png_files, 1):
            content_parts.append(f"![Figure {i}]({os.path.join(args.input, 'files', png_file)})\n*Figure {i}: {os.path.splitext(png_file)[0]}*\n\n")

    # Add Scripts Analysis section
    if scripts_files:
        content_parts.append("## Scripts Analysis\n\n")
        for script_file in scripts_files:
            with open(os.path.join(args.input, 'scripts', script_file), 'r') as f:
                script_content = f.read()
            
            # Call DeepSeek to analyze script
            script_payload = {
                "model": "deepseek-chat",
                "messages": [
                    {
                        "role": "system",
                        "content": "Analyze this script and provide a concise method description for all the analysis if applicable."
                    },
                    {
                        "role": "user", 
                        "content": f"Script: {script_file}\n\n{script_content}"
                    }
                ],
                "temperature": 0.5
            }
            
            script_response = requests.post(DEEPSEEK_API_URL, json=script_payload, headers=headers)
            script_response.raise_for_status()
            
            content_parts.append(f"### {script_file}\n")
            content_parts.append(script_response.json()["choices"][0]["message"]["content"] + "\n\n")

    # Combine all content
    full_content = "# Project Documentation\n\n" + "".join(content_parts)
    full_content += f"Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    
    # Call DeepSeek API to generate final markdown
    headers = {
        "Authorization": f"Bearer {DEEPSEEK_API_KEY}",
        "Content-Type": "application/json"
    }
    
    payload = {
        "model": "deepseek-chat",
        "messages": [
            {
                "role": "system",
                "content": "You are a helpful assistant that generates well-formatted markdown documentation."
            },
            {
                "role": "user",
                "content": f"Please format this content as a bioinformatics markdown in Chinese. Fig and caption should be separated by an empty line. Summarise analysis method in the last section as a bullet list:\n\n{full_content}"
            }
        ],
        "temperature": 0.7
    }
    
    response = requests.post(DEEPSEEK_API_URL, json=payload, headers=headers)
    response.raise_for_status()
    
    # Save the generated markdown
    with open(args.output+'.md', "w") as f:
        f.write(response.json()["choices"][0]["message"]["content"])
    
    print(f"Markdown documentation generated: {args.output}")
    
if __name__ == "__main__":
    generate_markdown_with_deepseek()