![png](https://github.com/yemingx/xieyeming1.github.io/blob/master/blood1.png)

| point | note | recommend |
| --- | --- | --- |
| extraction | gDNA/cfDNA/whole blood | cfDNA + WBC gDNA |
| seq | WGBS/target | target BS |
| QC | WBC/ctDNA | WBC filter + ctDNA marker |
| features | cpg/snp/nucleosome | CpG+SNP+nucleosome |
| cohort_EDA | disease/age/lifestyle/sex | AI agent recording |
| analysis | Elastic Net Regression | multimodal |
| benchmark | cor/MAE | age cor + MAE + organ age |
| application | aging | aging+organ_aging+MRD |


Boarder scope:
- 跨部门协作沟通，让colleague理解项目本身的价值和意义，任务拆解，找到卡点安排优先级
- 长期数据管理，生信标准化平台
- LDT（携带者hic，poreC）
- 整合到抗衰老等药物研发管线
- 药物新靶点预测



### 做cfDNA甲基化和SNP分析的时候，什么场景下，比如MRD，要考虑造血克隆污染？如何去除噪音，一定需要同时测一管血离心后的血浆cfDNA和白细胞吗?

#### 在基于血浆cfDNA的甲基化和SNP分析中，克隆性造血（Clonal Hematopoiesis, CH）的干扰是影响检测准确性的关键问题，尤其在微小残留病灶（MRD）监测、癌症早筛和低频突变分析等场景中需重点考虑。是否需要同步检测白细胞（WBC）取决于具体应用场景和技术策略，以下结合最新研究进展进行系统分析：
---
### 一、需考虑克隆造血（CH）干扰的核心场景
#### 1. MRD（微小残留病灶）监测
- 干扰机制：
MRD检测需识别极低丰度的肿瘤信号（ctDNA VAF常<0.1%），而CH相关突变（如DNMT3A、TET2）在cfDNA中的VAF通常≥2%。若无WBC对照，CH突变易被误判为肿瘤残留信号，导致假阳性。
- 临床影响：
在结直肠癌和乳腺癌的MRD研究中，未过滤CH的cfDNA检测假阳性率可升高50%以上，错误提示复发风险。
#### 2. 癌症早筛与多癌种检测
- 干扰机制：
早筛依赖低频突变信号（如TP53、KRAS），但CHIP携带者中约10%存在相同基因突变（如TP53在CH中占5-10%）。
- 案例：
CancerSEEK技术因未整合WBC过滤，早期癌症灵敏度仅43%。
#### 3. 低频SNP或突变分析
- 干扰机制：
肿瘤驱动基因（如JAK2、ASXL1）在CH和肿瘤中均高频突变，且CH突变在cfDNA中占比高达50%。
---
### 二、去除CH干扰的方法与技术策略
#### 1. 金标准：同步检测血浆cfDNA和WBC
- 原理：
通过配对WBC测序，识别并过滤cfDNA中与WBC共享的突变（即CH来源）。
- 操作流程：
    - 同管血离心分离血浆（cfDNA源）和Buffy Coat（WBC源）；
    - WBC与cfDNA同步进行高深度测序（建议≥500×）；
    - 生物信息学过滤：保留仅存在于cfDNA且WBC中VAF<0.1%的突变。
- 优势：
可降低73%的假阳性率，尤其适用于MRD和早筛。
#### 2. 无需WBC的替代方案
> 注：此类方法适用于无法获取WBC的场景，但精度可能受限。
- 甲基化特征分析：
    - CH不影响甲基化模式，因甲基化具有组织特异性（如肝癌特异性标志物RASSF1A）。
    - 应用：结直肠癌MRD检测中，ctDNA甲基化标志物（WIF1/NPY）的敏感性与突变联合分析相当，且无需WBC对照。
- 片段组学筛选：
    - 肿瘤ctDNA片段更短（中位长度144 bp vs. CH来源cfDNA的167 bp）；
    - 通过片段大小分布模型区分肿瘤与CH信号。
- 机器学习模型：
    - 基于突变特征（如CH偏好DNMT3A/TET2）、片段末端序列等训练分类器；
    - 例如：随机森林模型仅用cfDNA数据区分CH与肿瘤突变（AUC=0.90）。
- VAF阈值与基因背景过滤：
    - 排除VAF>2%的突变（因肿瘤早期VAF<1%）；
    - 过滤非驱动基因或非肿瘤相关突变（如CH高频基因DNMT3A）。
---
### 三、不同场景下的优化策略与临床建议
|应用场景|是否必需WBC|推荐技术组合|临床证据等级|
| --- | --- | --- | --- |
|MRD监测|强烈推荐|WBC过滤 + ctDNA甲基化|1级|
|癌症早筛|推荐|WBC过滤 + 片段组学 + VAF阈值|2级|
|治疗反应动态监测|可选|ctDNA甲基化（无需WBC）|2级|
|肿瘤溯源|非必需|组织特异性甲基化谱|3级|

---
### 四、总结：WBC检测的必要性与替代方案
1. 必需WBC的场景：
    - MRD监测、低频突变早筛、突变丰度低（VAF<1%）的分析；
    - 经济允许时优先选择（成本较高但假阳性率最低）。
2. 无需WBC的方案：
    - 甲基化分析：适用于治疗反应监测、肿瘤溯源；
    - 片段组学/机器学习：适用于早筛初筛或资源有限场景；
    - VAF过滤：作为辅助手段，但可能遗漏低频肿瘤信号。
> 临床建议：根据美国国立卫生研究院（NIH）指南，MRD相关研究需强制报告WBC对照使用情况；早筛产品（如Grail的Galleri）已逐步整合多组学特征以减少对WBC的依赖。

