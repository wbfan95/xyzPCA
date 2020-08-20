# xyzPCA
A Fortran code of principal component analysis (PCA). 

- **Author** : Wenbin, FAN (fanwenbin@shu.edu.cn)
- **Supervisor** : Yongle, LI

# Usage
## Prerequisite
1. Intel Fortran Compiler, OpenMP are needed. 
2. Only `.xyz` format is supported. Other format could be converted using VMD. 
3. `Make` only once. Then `./xyzPCA <Your .xyz file>`. 
## Output
1. Average (mean) structure in `*_mean.xyz` file. 
2. All PC values with percentages in `*_PCvalue` file. 
3. First five PC vectors in `*_PCvector` file, cooresponding to first five PC values. 
4. Fisrt three PCA in `*_PCA` file. 
