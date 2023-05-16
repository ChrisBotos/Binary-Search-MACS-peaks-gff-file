# file='C:\\Users\\Chris Botos\\PycharmProjects\\pythonProject\\Galaxy1-[Peaxi162annotation_v4_filtered.gff].gff'
# f=open(file)
# x=filesizer(file)
# t=peakfinder(f,'Peaxi162Scf00883',973127.5,0,3,4,'\t',x,5,5000,3000)

import time

start_time = time.time()

# Copy/paste the code into a python terminal and run it
def finder(file_name,MACS2,file_position_turn,file_start_turn,file_end_turn,divider,M_position_turn,M_start_turn,M_end_turn,dividerM,dig,b1,b2):
    r = []
    a = file_position_turn
    c = file_start_turn
    d = file_end_turn
    aM = M_position_turn
    cM = M_start_turn
    dM = M_end_turn
    # b1:bracket away from the gene(5000) b2:bracket towards the gene(3000) dig:number of digits at the end of chromosome name
    import os
    def filesizer(file_path):
        size = os.path.getsize(file_path)
        return (size)
    x = filesizer(file_name)
    def extender(line, b1, b2, c, d):
        if line[6] == '+':
            n = int(line[c]) - b1  # 5000
            y = int(line[c]) + b2  # 3000
        else:
            n = int(line[d]) - b2
            y = int(line[d]) + b1
        return (n, y)
    def peakfinder(f, p, peak, a, c, d, divider, x, dig, b1, b2):
        x = x // 2
        f.seek(x)
        f.readline()
        line = f.readline().strip().split(divider)
        tp = int(p[::-1][:dig][::-1])
        m = 0
        beginning = True
        k = int(x)
        t = None
        failsafe = 0
        while (p != line[a]):  # find the position/chromosome
            t = int(line[a][::-1][:dig][::-1])
            if failsafe > 13:  # failsafe
                h = 0
                if t > tp:
                    if not beginning:
                        f.seek(m)
                        f.readline()
                        while (p != line[a]):
                            x = f.tell()
                            line = f.readline().strip().split(divider)
                    else:
                        f.seek(0)
                        while (p != line[a]):
                            x = f.tell()
                            line = f.readline().strip().split(divider)
                else:
                    while (p != line[a]):
                        x = f.tell()
                        line = f.readline.strip().split(divider)
            if t > tp:
                k = k // 2
                x = x - k
                f.seek(x)
                f.readline()
                line = f.readline().strip().split(divider)
            else:
                beginning = False
                m = x
                k = k // 2
                x = x + k
                f.seek(x)
                f.readline()
                line = f.readline().strip().split(divider)
            failsafe = failsafe + 1
        failsafe2 = 0
        n = extender(line, b1, b2, c, d)[0]
        y = extender(line, b1, b2, c, d)[1]
        while not (n <= peak and peak < y):  # find the exact position of the peak
            if failsafe2 > 13:  # failsafe
                h = 0
                if peak > y:
                    if not beginning:
                        f.seek(m)
                        f.readline()
                        line = f.readline().strip().split(divider)
                        n = extender(line, b1, b2, c, d)[0]
                        y = extender(line, b1, b2, c, d)[1]
                        while not (n <= peak and peak < y) or p != line[a]:
                            x = f.tell()
                            line = f.readline().strip().split(divider)
                            n = extender(line, b1, b2, c, d)[0]
                            y = extender(line, b1, b2, c, d)[1]
                            h = h + 1
                            if h > k:
                                return (x, ['Not aligned'])
                    else:
                        f.seek(0)
                        while not (n <= peak and peak < y) or p != line[a]:
                            x = f.tell()
                            line = f.readline().strip().split(divider)
                            n = extender(line, b1, b2, c, d)[0]
                            y = extender(line, b1, b2, c, d)[1]
                            h = h + 1
                            if h > k:
                                return (x, ['Not aligned'])
                else:
                    while not (n <= peak and peak < y) or p != line[a]:
                        x = f.tell()
                        line = f.readline().strip().split(divider)
                        n = extender(line, b1, b2, c, d)[0]
                        y = extender(line, b1, b2, c, d)[1]
                        h = h + 1
                        if h > k:
                            return (x, ['Not aligned'])
                return (x, line)
            if peak < y:
                k = k // 2
                x = x - k
                f.seek(x)
                f.readline()
                line = f.readline().strip().split(divider)
                n = extender(line, b1, b2, c, d)[0]
                y = extender(line, b1, b2, c, d)[1]
                while (line[a] != p):
                    k = k // 2
                    x = x + k
                    f.seek(x)
                    f.readline()
                    line = f.readline().strip().split(divider)
                    n = extender(line, b1, b2, c, d)[0]
                    y = extender(line, b1, b2, c, d)[1]
                    if k <= 1:
                        failsafe2 = 23
                        break
            if n <= peak:
                beginning = False
                m = x
                k = k // 2
                x = x + k
                f.seek(x)
                f.readline()
                line = f.readline().strip().split(divider)
                n = extender(line, b1, b2, c, d)[0]
                y = extender(line, b1, b2, c, d)[1]
                while (line[a] != p):
                    print('d', k)
                    k = k // 2
                    x = x - k
                    m = x
                    f.seek(x)
                    f.readline()
                    line = f.readline().strip().split(divider)
                    if k <= 1:
                        failsafe2 = 23
                        break
            failsafe2 = failsafe2 + 1
        return (x, line)  # returns the byte position the line after which contains the chosen gene ID and that line
    with (open(file_name) as f, open(MACS2) as g):
        h=0
        for line in g:
            line=line.strip().split(dividerM)
            p=line[aM]
            peak=float(line[cM])+(float(line[dM])-float(line[cM]))/2
            t=peakfinder(f,p,peak,a,c,d,divider,x,dig,b1,b2)
            r=r+[t[1]]
            h=h+1
            print(h)
    def cleanser(r):
        cl = []
        while ['Not aligned'] in r:
            r.remove(['Not aligned'])
        for i in r:
            if i[8][:3] == 'ID=':
                cl = cl + [i[8][3:27]]
            else:
                cl = cl + [i[8][7:31]]
        return (cl)
    def peaks_per_gene(cl):  # Score behind gene
        l = []
        for i in cl:
            if i not in l:
                l = l + [1]
                l = l + [i]
            else:
                l[l.index(i) - 1] = l[l.index(i) - 1] + 1
        return (l)
    res=peaks_per_gene(cleanser(r))
    for i in res:
        if type(i)==int:
            continue
        if i[::-1][:2]=='N;':
            res[res.index(i)]=i[:len(i)-2]
    return(res)



# o=finder('Galaxy1-[Peaxi162annotation_v4_filtered.gff].gff','MACS2Ypol.txt',0,3,4,'\t',0,1,2,'\t',5,5000,3000)
# print(o)
# The number in front of the gene shows the number of peaks found in it

end_time = time.time()
total_time = end_time - start_time
print("Total time taken: ", total_time, "seconds")

