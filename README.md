# 1.LLM application for bioinformatics workflow

### Achieved
#### email_fetch (测序下机邮件信息抓取)
https://github.com/yemingx/xieyeming1.github.io/tree/master/llm_app/email_fetch
#### wechatPub_fetch （微信公众号文章信息抓取）
https://github.com/yemingx/xieyeming1.github.io/tree/master/llm_app/wechatPub

### Ongoing
#### agent, 输入markdown，输出method和fig legend
#### epigenome/3d genome/transcriptome agent
表观数据标准，输入提示词和清洗好的表观，输出代码，表格和图（png pdf svg）；提供一个教程，有详细的讲解和实例，做一个数据进行可视化EDA的agent，使用graphQL存储和调用已有的生物信息学多组学表格。现有的处理好的数据为表观基因组数据的peak和有信号值的bedgraph文件，包括dnase seq，atac seq，各种转录因子的chip seq。agent可以通过提示词输入，画出peak overlap的eulerr plot, upset plot，能在图上标注必备的统计数值，另外也能调用本地安装的deeptools画peak信号的堆积图。
- agent, 输入narrowpeak；跑 R；得到upset plot, eulerr plot

# 2.GPU acceleration 
### TBD
#### Parabricks for bowtie2 bwa 
#### hic analysis
#### C++基础→SIMD优化（SSE/AVX指令集）→CUDA并行计算基础; 《算法导论》核心章节→实现经典生信算法（如BLAST优化、序列比对算法复现）; CMake构建系统 → GitHub CI/CD流程 → 打包Python/R工具链。 CUDA原子操作/共享内存 → cuDNN内核定制 → PyTorch算子开发

# 3.AIDD
### TBD
#### learning ML/DL
#### AI药筛; AI设计binder

# 4.AI modeling project
### TBD
#### gemma3 SLM fine-tune; google/gemma-3-270m
#### atac（dnase）to MICC diffusion 文生图
#### dipC 12878 to hicirc 12878 diffusion
