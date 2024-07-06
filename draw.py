import os
from numpy import *
import matplotlib.pyplot as plt
import sys

def show_and_save(rho, output_file):
    x = linspace(-0.1, 3.1, 301)
    y = linspace(-0.1, 1.1, 101)
    X, Y = meshgrid(x, y)
    plt.figure()
    plt.contourf(X, Y, rho)
    plt.colorbar()
    plt.title('t=4.0')
    plt.xlabel("x")
    plt.savefig(output_file)  # 保存图片
    plt.close()

# 修改主程序部分，循环处理多个文件
for i in range(1, 6):  # 假设要处理output1到output10这十个文件
    filename = f"output{i}.txt"
    if os.path.exists(filename):
        rho = loadtxt(filename)
        # 这里可以添加一些数据处理，如你之前注释掉的代码

        # 调用show_and_save函数显示并保存图像
        output_image = f"output{i}_plot.png"
        show_and_save(rho, output_image)
    else:
        print(f"文件 {filename} 不存在。")

print("所有文件处理完成。")
