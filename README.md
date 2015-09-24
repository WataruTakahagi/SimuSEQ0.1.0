## <span style="color:red">S</span>imuSE<span style="color:blue">Q</span>

### Install
```shell
python setup.py install
# or user install with
# python setup.py install --user
```

### Microarray analysis tool for whole-ecoli project

    from simuseq import SimuSEQ  
    SimuSEQ().auto([1,1,1,1,0,1])
    
    > python main.py
    SimuSEQ().readseq()    : on  
    SimuSEQ().makeprobe()  : on  
    SimuSEQ().microarray() : on  
    SimuSEQ().rank()       : on  
    SimuSEQ().fasta()      : off  
    SimuSEQ().makedata()   : on  
    Continue? (y/n) : y
    
![microarray.png](http://www.fastpic.jp/images.php?file=8543443149.png "microarray.png")
    
If you would like to check the status of SimuSEQ, please run main.py .  

    > python main.py
    SimuSEQ().readseq()    : on  
    SimuSEQ().makeprobe()  : on  
    SimuSEQ().microarray() : on  
    SimuSEQ().rank()       : on  
    SimuSEQ().fasta()      : off  
    SimuSEQ().makedata()   : on  
    Continue? (y/n) : n
    
You can change arguments to set the status of SimuSEQ auto mode
    
### 1 : on, 0 : off
#### SimuSEQ( ).auto([1,1,1,1,0,1])
1: SimuSEQ( ).readseq( )  
1: SimuSEQ( ).makeprobe( )    
1: SimuSEQ( ).microarray( )    
1: SimuSEQ( ).rank( )        
0: SimuSEQ( ).fasta( )  
1: SimuSEQ( ).makedata( )    

