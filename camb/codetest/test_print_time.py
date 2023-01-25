import datetime
 
# 打印当前时间
time1 = datetime.datetime.now()
print(time1)
# 打印按指定格式排版的时间
time2 = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print(time2)
