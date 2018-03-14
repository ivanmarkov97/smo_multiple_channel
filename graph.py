import matplotlib.pyplot as plt
import re

data_files = [
        "file_oper_refuse.txt",
        "file_load_koef.txt",
        "file_refuse_reverse.txt",
        "file_load_koef2.txt",
        "file_len_queue2.txt",
        "file_work_queue2.txt",
        "file_time_queue2.txt",
        "file_load_koef4.txt",
        "file_len_queue.txt",
        "file_time_queue.txt",
        "file_load_koef5.txt",
        "file_len_queue5.txt",
        "file_time_queue5.txt",
        "file_free_works.txt",
        "file_work_koef.txt"        
]

for df in data_files:
        m_file = open(df)
        X = []
        Y = []
        for r in m_file:
                res = re.findall(r"[-+]?\d*\.\d+|\d+", r)
                print res
                X.append(float(res[0]))
                Y.append(float(res[1]))
        plt.plot(X, Y)
        plt.ylabel(df)
        plt.show()
