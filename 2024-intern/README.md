# 2024-intern-exam

## Instructions

This is an open-book assessment that may take a couple hours to complete. All tools, languages, and resources are fair game. **Fork this repository.** The exam submission should be made as a commit to your forked repo.

## Prompts

#### 1) Debug this code
The code in this directory performs the core preprocessing workflow of the scanpy framework, but it does not currently work.
   
   - Assess the code for issues, and contribute changes that will allow it to run to completion according to single cell processing best practices.
   - Elaborate on your changes and logic with inline comments.

A conda `environment.yml` file that accounts for necessary dependencies is included for convenience. Only code changes submitted to your forked repo are required. **Don't** submit the anndata object.

```
# when refactored correctly, the workflow should execute with this command
python main.py --output_file "out.h5ad"
```

#### 2) Consider the unknown sequencing data given [here](https://drive.google.com/drive/folders/15AL1nuJCV2EC9p0LIl0YudXVraAjsOay?usp=drive_link).
   - Perform a QC process on this data.
     - a suggested approach would be to use a tool combo like fastQC/MultiQC
   - Provide a detailed assessment of the nature and quality of this data given your QC output.
   - Speculate on the type of sequencing used to generate this data.
   - Summarize your conclusions as a markdown file in your forked repo. If you have relevant outputs like a MultiQC summary, please include them.
