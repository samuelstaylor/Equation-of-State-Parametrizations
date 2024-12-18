import os

def count_lines_in_repo(directory, extension=".py"):
    total_lines = 0
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith(extension):
                file_path = os.path.join(dirpath, filename)
                with open(file_path, 'r', encoding='utf-8') as f:
                    total_lines += sum(1 for line in f if line.strip())  # Count non-blank lines
    return total_lines

if __name__ == "__main__":
    repo_path = "."  # Current directory (your repository root)
    lines_of_code = count_lines_in_repo(repo_path)
    print(f"Total number of lines in Python files: {lines_of_code}")
