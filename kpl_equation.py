
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 15:08:54 2018

@author: ysc
"""

#导入各种模块



import numpy as np  #导入numpy模块,并命名为np
import math
import matplotlib #导入matplotlib模块用于绘图
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab! #在jupyter notebook上运行要注释掉此句 

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from moviepy.video.io.bindings import mplfig_to_npimage
import moviepy.editor as mpy


import sympy as sp #导入sympy包中的函数


from mpl_toolkits.mplot3d import Axes3D

from scipy import integrate





e_ = 0.6#离心率
p = 6 #半通经

r1 = p/(1-e_)#远心距
r2 = p/(1+e_)#近心距


a_ = (1.0/2.0) * (r1 + r2)
c_ = a_ - r2
b_ = math.sqrt(a_**2 - c_**2)


omega_r = 5 #平均角速度
omega_phi = 0.0#平均角速度

T_r = 2 * math.pi / omega_r #回归周期

M_0 = 0#平近点角初相
M2_0 = 0#进动平近点角初相

def E_M_(M_):
    if M_==0:
        return 0
    elif M_ == 2 * math.pi:
        return 2 * math.pi
    else:
        max_ = 2 * math.pi  
        min_ = 0.0

        count=1 
        while (max_-min_)>1e-8:
            med_ = (max_ + min_) / 2.0   
            a = min_ - e_ * math.sin(min_) - M_
            b = max_ - e_ * math.sin(max_) - M_
            c = med_ - e_ * math.sin(med_) - M_
            count+=1
            if  a * c <0: 
                max_ = (min_ + max_)/2.0
            else:
                if b * c <0:
                    min_ = (min_ + max_)/2.0
                else:
                    break
        return ((min_ + max_)/2.0)


def multi(M_):
    if(math.sin(M_)>=0):
        return math.acos(math.cos(M_))
    else:
        return 2*math.pi - math.acos(math.cos(M_))
       

#绘图$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# ta = 0 #初始时刻 
# tb = 2 * math.pi#终止时刻
# step = 0.1 #步长
# number = int(tb / step)

# x_plot = [0 for x_ in range(number)]
# y_plot = [0 for x_ in range(number)]


# for i in range(1, number):
#     x_plot[i] = x_plot[i-1] + step
#     y_plot[i] = E_M_(x_plot[i])


# #二维xy绘图 
# plt.figure(figsize = (6,6))
# plt.plot(x_plot[1:number], y_plot[1:number])
# # plt.plot(x_plot, y2_plot)
# plt.ylabel('Y')#坐标轴
# plt.xlabel('X')
# # plt.xlim(0,2*math.pi)
# # plt.ylim(0,2*math.pi)

# plt.show()#显示图片
# plt.savefig('graph.png',format = 'png', dpi = 300)
# plt.close()#关闭图片


# #绘动图$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


phi_p = np.arange(0,2 * math.pi,0.01)#创建列表以存储坐标
r_p = p / (1 + e_ * np.cos(phi_p))
x_p = r_p * np.cos(phi_p)
y_p = r_p * np.sin(phi_p)


# 使用MATPLOTLIB画一个图
duration = 15#总共的时间       figsize代表横纵比
fig_mpl, ax = plt.subplots(1, figsize=(6,6),facecolor='white')
r1 = p/(1-e_) 
r2 = p/(1+e_)
x1 = r1 * math.cos(0)
y1 = r1 * math.sin(0)

x2 = r2 * math.cos(math.pi)
y2 = r2 * math.sin(math.pi)






ax.set_xlim(-20,20)#设置纵坐标
ax.set_ylim(-20,20)#设置纵坐标


line = ax.plot(x1,y1,'o',ms=1, c ='black')#远点
line = ax.plot(x1,y1,'o',ms=1, c ='black')#远点
line = ax.plot(x1,y1,'o',ms=1, c ='black')#远点
# line = ax.plot(x1,y1,'o',ms=1, c ='black')#远点
# line = ax.plot(x1,y1,'o',ms=1, c ='black')#远点
# line = ax.plot(x1,y1,'o',ms=1, c ='black')#远点
# line = ax.plot(x1,y1,'o',ms=1, c ='black')#远点

# 使用MOVIEPY让图动起来(根据时间t来更新图). 保存为GIF.
def make_frame_mpl(t):
    n = 8 #几圈
    t_ = n * T_r * t / duration


    M_ = omega_r * t_ + M_0#平近点角
    M2_ = omega_phi * t_ + M2_0#进动平近点角


    
    M_m = multi(M_)
    M2_m = - multi(M2_)
    # M_m = math.asin(math.sin(M_))#多个周期

    #粒子轨迹
    xi_ =  E_M_(M_m)#偏近点角

    x_plota = a_ * (math.cos(xi_) - e_)
    y_plota = a_ * (math.sqrt(1 - e_**2)) * math.sin(xi_) 

    x_plot2a = x_plota
    y_plot2a = y_plota / math.sqrt(1 - e_**2)


    ax.lines.pop()#不留轨迹
    ax.lines.pop()#不留轨迹
    ax.lines.pop()#不留轨迹
    # ax.lines.pop()#不留轨迹
    # ax.lines.pop()#不留轨迹
    # ax.lines.pop()#不留轨迹
    # ax.lines.pop()#不留轨迹


    phi_p = np.arange(0,2 * math.pi,0.01)#创建列表以存储坐标
    r_p = p / (1 + e_ * np.cos(phi_p)) #椭圆方程


    x_pa = r_p * np.cos(phi_p) #椭圆
    y_pa = r_p * np.sin(phi_p)

    x_p2a = x_pa #辅助圆
    y_p2a= y_pa / math.sqrt(1 - e_**2)

#旋转
    x_plot = x_plota * math.cos(M2_m) + y_plota * math.sin(M2_m)
    y_plot = y_plota * math.cos(M2_m) - x_plota * math.sin(M2_m)

    x_plot2 = x_plot2a * math.cos(M2_m) + y_plot2a * math.sin(M2_m)
    y_plot2 = y_plot2a * math.cos(M2_m) - x_plot2a * math.sin(M2_m)


    x_p = x_pa * math.cos(M2_m) + y_pa * math.sin(M2_m)
    y_p = y_pa * math.cos(M2_m) - x_pa * math.sin(M2_m)

    x_p2 = x_p2a * math.cos(M2_m) + y_p2a * math.sin(M2_m)
    y_p2 = y_p2a * math.cos(M2_m) - x_p2a * math.sin(M2_m)

 
    # line = ax.quiver(x_plot, y_plot)#粒子轨迹方向


    line = ax.plot(x_plot, y_plot,'o',ms=1.0 , c ='black')#粒子轨迹 
    # line = ax.plot(x_plot2, y_plot2,'o',ms=1.0 , c ='black')#像粒子轨迹
    line = ax.plot(x_plot, y_plot,'o',ms=5 , c ='red')#粒子
    # line = ax.plot(x_plot2, y_plot2,'o',ms=5 , c ='blue')#像粒子
    line = ax.plot(x_p, y_p,'o',ms=0.1 , c ='red')#椭圆
    # line = ax.plot(x_p2, y_p2,'o',ms=0.1 , c ='blue')#辅助圆
    line = ax.plot(0,0,'o',ms=1, c ='red')#原点
    # line = ax.plot(-c_,0,'o',ms=1, c ='blue')#椭圆中心
    # line = ax.plot(r1 * math.cos(0+2*math.pi *f * t1/duration),r1 * math.sin(0+2*math.pi *f * t1/duration),'o',ms=3, c ='black')#远点

    ax.set_title("t = " + str(round(T_r * t / duration,1)))
    return mplfig_to_npimage(fig_mpl) # RGB image of the figure

animation =mpy.VideoClip(make_frame_mpl, duration=duration)
animation.write_gif("tests2D.gif", fps=30)


