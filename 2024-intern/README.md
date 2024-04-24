# 2024-intern-exam

## Instructions

- This is an open-book assessment that may take a couple hours to complete.
- All tools, languages, and resources are fair game.
- **Fork this repository.** The exam submission should be made as a commit to your forked repo.

## Prompts

1. **Debug this code**
   - The code in this directory perform the core preprocessing workflow of the scanpy framework, but it does not currently work.
   - Assess the code for issues and contribute changes that will allow it to run to completion according to single cell preprocessing best practices.
   - A conda environment.yml file is included for convenience.
1. **Consider the unknown sequencing data given [here](https://drive.google.com/drive/folders/15AL1nuJCV2EC9p0LIl0YudXVraAjsOay?usp=drive_link).**
   - Perform a QC process on this data, and aggregate your results into a single output.
     - a suggested approach would be to use a tool combo like fastQC/MultiQC
   - Provide a detailed assessment of the nature and quality of this data given your QC output.
   - Speculate on the type of sequencing used to generate this data.
