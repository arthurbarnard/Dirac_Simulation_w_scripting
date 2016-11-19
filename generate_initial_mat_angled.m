dX=2.5;
Nx=1201;
Ny=721;
Xoffset=1500;
Yoffset=0;
[X,Y]=meshgrid((-(Nx-1)):(Nx),(-(Ny-1)):(Ny));
X=X*dX/2+Xoffset;
Y=Y*dX/2+Yoffset;

Vx=mysmoothbarrier(X(1,:),1390,20,1560,20);

V=ones(Ny*2,1)*Vx;


AbsMat=zeros(size(X));
AbsMat(abs(Y)>800)=1;
AbsMat(X>2500)=1;

DriveMat=zeros(size(X));
DriveMat((X>-5)&(X<10)&(Y>-150)&(Y<150))=1;

save('Klein_barrier_20nm_edges.mat','dX','Nx','Ny','Xoffset','Yoffset','V','AbsMat','DriveMat');

