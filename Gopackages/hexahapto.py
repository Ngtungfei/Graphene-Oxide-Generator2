import numpy as np
import matplotlib.pyplot as plt

def buildhexahapto(n):
    step=1
    cabond=1.34
    carbonlines=[]
    atomsxyza=[]
    atomsxyz=[]
    bondsxyz=[]
    atomno=0
    xmin=0
    ymin=0
    xmax=0
    ymax=0
    
    for i in range(int(n//2+1)):
        atomsxyza,carbonlines=hax(i*2+1,atomsxyza,carbonlines)
           
    for i in range(len(carbonlines)):
        atomsxyz.append(str(i+1))
        atomsxyz.append('C')
        atomsxyz.append(atomsxyza[i*2])
        atomsxyz.append(atomsxyza[i*2+1])
        atomsxyz.append(0)
        atomsxyz.append('CG2R61')
        atomsxyz.append('1')
        atomsxyz.append('GO')
        atomsxyz.append('0.0')
    aa=0
    lines=[]
    for i in range(int(n//2+1)):
        a=int(6*((i-1)*2+1))
        if a<0:
            a=0
        aa+=a
        lines.append([])
        for j in range(int(6*(i*2+1))):
            lines[i].append(j+1+aa)
            if j+1==int(6*(i*2+1)):
                bondsxyz.append(int(6*(i*2+1))-j+aa)
                bondsxyz.append(j+1+aa)
                bondsxyz.append(int((i+j)%2+1))
            else:
                bondsxyz.append(j+1+aa)
                bondsxyz.append(j+2+aa)
                bondsxyz.append(int((i+j)%2+1))
                 
    npatomsdata=np.array(atomsxyz).reshape(int(len(atomsxyz)/9),9)
    npbondsdata=np.array(bondsxyz).reshape(int(len(bondsxyz)/3),3)
    for i in range(int(n//2+1)-1):
        for a in lines[i]:
            for b in lines[i+1]:
                x=(float(npatomsdata[a-1][2])-float(npatomsdata[b-1][2]))**2+(float(npatomsdata[a-1][3])-float(npatomsdata[b-1][3]))**2
                if x < 1.8:
                    bondsxyz.append(a)
                    bondsxyz.append(b)
                    bondsxyz.append(1)
       
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

def hax(n,atomsxyz,carbonlines):
    kk=0
        
    cabond=1.34
    atomsxyz.append(-n*1.732*cabond/2)
    atomsxyz.append(cabond/2)
    kk+=1
    carbonlines.append(kk)
        
    for i in range(int(n//2)):
        atomsxyz.append(atomsxyz[-2]+0.866*cabond)
        atomsxyz.append(atomsxyz[-2]+0.5*cabond)
        kk+=1
        carbonlines.append(kk)
        atomsxyz.append(atomsxyz[-2])
        atomsxyz.append(atomsxyz[-2]+cabond)
        kk+=1
        carbonlines.append(kk)
       
    for i in range(int((n//2+1)*2-1)):
        if i%2==0:
            atomsxyz.append(atomsxyz[-2]+0.866*cabond)
            atomsxyz.append(atomsxyz[-2]+0.5*cabond)
            kk+=1
            carbonlines.append(kk)
        else:
            atomsxyz.append(atomsxyz[-2]+0.866*cabond)
            atomsxyz.append(atomsxyz[-2]-0.5*cabond)
            kk+=1
            carbonlines.append(kk)
            
    for i in range(int(n//2)):
        atomsxyz.append(atomsxyz[-2]+0.866*cabond)
        atomsxyz.append(atomsxyz[-2]-0.5*cabond)
        kk+=1
        carbonlines.append(kk)
        atomsxyz.append(atomsxyz[-2])
        atomsxyz.append(atomsxyz[-2]-cabond)
        kk+=1
        carbonlines.append(kk)

    atomsxyz.append(n*1.732*cabond/2)
    atomsxyz.append(cabond/2)
    kk+=1
    carbonlines.append(kk)
    atomsxyz.append(n*1.732*cabond/2)
    atomsxyz.append(-cabond/2)
    kk+=1
    carbonlines.append(kk)

    for i in range(int(n//2)):
        atomsxyz.append(atomsxyz[-2]-0.866*cabond)
        atomsxyz.append(atomsxyz[-2]-0.5*cabond)
        kk+=1
        carbonlines.append(kk)
        atomsxyz.append(atomsxyz[-2])
        atomsxyz.append(atomsxyz[-2]-cabond)
        kk+=1
        carbonlines.append(kk)

    for i in range(int((n//2+1)*2-1)):
        if i%2==0:
            atomsxyz.append(atomsxyz[-2]-0.866*cabond)
            atomsxyz.append(atomsxyz[-2]-0.5*cabond)
            kk+=1
            carbonlines.append(kk)
        else:
            atomsxyz.append(atomsxyz[-2]-0.866*cabond)
            atomsxyz.append(atomsxyz[-2]+0.5*cabond)
            kk+=1
            carbonlines.append(kk)

    for i in range(int(n//2)):
        atomsxyz.append(atomsxyz[-2]-0.866*cabond)
        atomsxyz.append(atomsxyz[-2]+0.5*cabond)
        kk+=1
        carbonlines.append(kk)
        atomsxyz.append(atomsxyz[-2])
        atomsxyz.append(atomsxyz[-2]+cabond)
        kk+=1
        carbonlines.append(kk)

    atomsxyz.append(-n*1.732*cabond/2)
    atomsxyz.append(-cabond/2)
    kk+=1
    carbonlines.append(kk)
        
    return atomsxyz,carbonlines
def main():
    dic={1:'k',2:'r',3:'b',4:'g'}
    col={'C':'r','O':'b','H':'g'}
    siz={'C':'600','O':'650','H':'300'}

    npatomsdata,npbondsdata,xmin,xmax,ymin,ymax=buildhexahapto(5)
    #plt.xlim([xmin-4,xmax+4])
    #plt.ylim([ymin-3,ymax+3])
    for i in npbondsdata:
        plt.plot([float(npatomsdata[i[0]-1][2]),float(npatomsdata[i[1]-1][2])],[float(npatomsdata[i[0]-1][3]),float(npatomsdata[i[1]-1][3])],color=dic[i[2]])
    plt.axes().set_aspect('equal')
    #plt.title('%d Carbons %d left:\n %d(=O) %d(COOH) %d(-O-) %d(OH)'%(self.cno,self.clo,self.cono,self.coohno,self.ono,self.ohno))
    plt.show()

if __name__ == '__main__':
    main()  

    
