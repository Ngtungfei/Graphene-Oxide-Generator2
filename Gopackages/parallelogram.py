import numpy as np
import matplotlib.pyplot as plt

def buiildarallelogram(n,m):
    cabond=1.34
    xmin=0
    ymin=0
    xmax=0
    ymax=0 
    carbonline0=[]
    carbonline1=[]
    atomsxyz=[]
    bondsxyz=[]
    
    for i in range(n):
        for j in range(m):
            atomsxyz.append(str(j+i*m+1))
            atomsxyz.append('C')
            carbonline0.append(j+i*m+1)
            if j%2==0:
                atomsxyz.append(j//2*3*cabond+i*cabond*1.5)
                atomsxyz.append(i*1.732*cabond/2)
                atomsxyz.append(0)
            else:
                atomsxyz.append(j//2*3*cabond+cabond+i*cabond*1.5)
                atomsxyz.append(i*1.732*cabond/2)
                atomsxyz.append(0)
            atomsxyz.append('CG2R61')
            atomsxyz.append('1')
            atomsxyz.append('GO')
            atomsxyz.append('0.0')
            
    npcarbonline0=np.array(carbonline0).reshape(int(len(carbonline0)/m),m)              
    for i in  range(npcarbonline0.shape[0]):
        if i==0:
            for j in range(int(m/2)):
                bondsxyz.append(npcarbonline0[i][j*2])
                bondsxyz.append(npcarbonline0[i][j*2+1])
                bondsxyz.append(2)
        else:
            for j in range(int(m/2)):
                bondsxyz.append(npcarbonline0[i][j*2])
                bondsxyz.append(npcarbonline0[i][j*2+1])
                bondsxyz.append(1)
                            
    for ii in  range(npcarbonline0.shape[0]-1):
        i=ii+1
        if i%2==0:
            for j in range(m-1):
                bondsxyz.append(npcarbonline0[i][j])
                bondsxyz.append(npcarbonline0[i-1][j+1])
                bondsxyz.append(2)
        else:
            for j in range(m-1):
                bondsxyz.append(npcarbonline0[i][j])
                bondsxyz.append(npcarbonline0[i-1][j+1])
                bondsxyz.append(1)
       
    npatomsdata=np.array(atomsxyz).reshape(int(len(atomsxyz)/9),9)
    npbondsdata=np.array(bondsxyz).reshape(int(len(bondsxyz)/3),3)
    
    for i in npatomsdata[:,2]:
        if xmin> float(i):
            xmin=float(i)
        if xmax< float(i):
            xmax=float(i)

    for i in npatomsdata[:,3]:
        if ymin> float(i):
            ymin=float(i)
        if ymax< float(i):
            ymax=float(i)

    return npatomsdata,npbondsdata,xmin,xmax,ymin,ymax

def main():
    dic={1:'k',2:'r',3:'b',4:'g'}
    col={'C':'r','O':'b','H':'g'}
    siz={'C':'600','O':'650','H':'300'}

    npatomsdata,npbondsdata,xmin,xmax,ymin,ymax=buiildarallelogram(29,16)
    plt.xlim([xmin-4,xmax+4])
    plt.ylim([ymin-3,ymax+3])
    for i in npbondsdata:
        plt.plot([float(npatomsdata[i[0]-1][2]),float(npatomsdata[i[1]-1][2])],[float(npatomsdata[i[0]-1][3]),float(npatomsdata[i[1]-1][3])],color=dic[i[2]])
    plt.axes().set_aspect('equal')
    #plt.title('%d Carbons %d left:\n %d(=O) %d(COOH) %d(-O-) %d(OH)'%(cno,clo,cono,coohno,ono,ohno))
    plt.show()

if __name__ == '__main__':
    main()  

    
