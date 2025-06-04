import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    # 그래프 시각화
    # 데이터 로드
    df = pd.read_csv("result/results_2025-06-04_16-29-11")

    # 출력 디렉토리 생성
    output_dir = "figures"
    os.makedirs(output_dir, exist_ok=True)

    # 시각화 스타일 설정
    sns.set(style="whitegrid")
    plt.rcParams["figure.figsize"] = (10, 6)

    # 정확도 변화
    def plot_accuracy_vs(variable):
        plt.figure()
        sns.lineplot(data=df, x=variable, y="BWT_accuracy", label="BWT", marker="o", color='orange')
        sns.lineplot(data=df, x=variable, y="brute_accuracy", label="Brute Force", marker="s", color='blue')
        plt.title(f"Accuracy vs {variable}")
        plt.xlabel(variable)
        plt.ylabel("Accuracy")
        plt.ylim(0, 1.0)
        plt.tight_layout()
        # 저장
        plt.savefig(f"{output_dir}/accuracy_vs_{variable}.png")
        plt.close()
    # 복원 시간 변화
    def plot_time_vs(variable):
        plt.figure()
        sns.lineplot(data=df, x=variable, y="BWT_time_sec", label="BWT", marker="o", color='orange')
        sns.lineplot(data=df, x=variable, y="brute_time_sec", label="Brute Force", marker="s", color='blue')
        plt.title(f"Reconstruction Time vs {variable}")
        plt.xlabel(variable)
        plt.ylabel("Time (seconds)")
        plt.tight_layout()
        # 저장
        plt.savefig(f"{output_dir}/time_vs_{variable}.png")
        plt.close()

    # 변수별로 시각화
    variables = ["genome_limit", "num_reads", "read_length", "mutation_rate", "max_mismatch"]
    for var in variables:
        plot_accuracy_vs(var)
        plot_time_vs(var)


if __name__ == '__main__':
    main()
