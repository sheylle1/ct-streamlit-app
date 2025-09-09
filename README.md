# Ct Value Extraction Tool

A Streamlit web application for extracting Ct values from PCR data CSV files.

## Features

- Processes CSV files from folders and subfolders
- Extracts Ct values for samples
- Provides positive control quality checks
- Fuzzy matching of sample names
- Export results to CSV or Excel format

## Installation

1. Clone this repository:
```bash
git clone https://github.com/your-username/ct-value-extraction-tool.git
cd ct-value-extraction-tool

Install dependencies:

bash
pip install -r requirements.txt
Run the application:

bash
streamlit run app.py
Usage
Enter the path to the folder containing your CSV files

Upload a sample list with an "Episode Number" column

Click "Extract Ct Values" to process the files

View results in the preview tabs

Download results in CSV or Excel format

Deployment on Streamlit Cloud
Fork this repository

Go to Streamlit Cloud

Click "New app" and connect your GitHub account

Select your repository and branch

Set the main file path to app.py

Click "Deploy"

text

4. **Initialize Git and push to GitHub**:

Open your terminal/command prompt and run:

```bash
# Navigate to your project folder
cd ct-value-extraction-tool

# Initialize git
git init

# Add all files
git add .

# Commit files
git commit -m "Initial commit: Ct Value Extraction Tool"

# Add your GitHub repository as remote
git remote add origin https://github.com/your-username/your-repo-name.git

# Push to GitHub
git branch -M main
git push -u origin main
Deploy to Streamlit Cloud:

Go to share.streamlit.io

Sign in with your GitHub account

Click "New app"

Select your repository, branch (main), and file path (app.py)

Click "Deploy"