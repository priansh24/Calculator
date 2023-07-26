#include<stdlib.h>
#include<stdio.h>
#include<math.h>

double result = 0;char exitprog = ' ';int counter=0;

double add()
{
    int n;  printf("Enter number of Inputs: "); scanf("%d",&n);
    int sum=0,a;    printf("\nEnter elements Separated by Enter key:\n");
    for(int i=1;i<=n;i++)
    {
        scanf("%d",&a);
        sum+=a;
    }
    return sum;
}
double subs()
{
    int n;  printf("Enter number of Inputs: "); scanf("%d",&n);
    double sum,sum2,a;    printf("\n\nEnter elements Separated by Enter key:\n"); scanf("%lf",&a);
    if(n==1)    return a;  
    else
    {
        sum=a;  sum2=a;
        for(int i=1;i<n;i++)
        {
            scanf("%lf",&a);
            sum-=a;
            sum2+=a;
        }
    }
    if(counter == 1)    return -sum;
    else                return sum2;
}
double multiply()
{
    if(counter == 1)    result = 1;
    int n;  printf("Enter number of Inputs: "); scanf("%d",&n);
    long int product=1;int a; printf("Enter elements Separated by Enter key:\n"); scanf("%d",&a);product=a;
    for(int i=1;i<n;i++)
    {
        scanf("%d",&a);
        product*=a;
    }
    return product;
}
double divide()
{
    double quo=0;int a;
    if(counter!=1)
    {
        quo = multiply();
        return quo;
    }
    result=1;
    int n;  printf("Enter number of Inputs: "); scanf("%d",&n);
    printf("Enter elements Separated by Enter key:\n"); scanf("%d",&a);  quo=a;
    if(n==1)    return (1/quo);
    else 
    {
        for(int i=1;i<n;i++)
        {
            scanf("%d",&a);
            quo/=a;
        }
    }
    return 1/quo;
}
double mod()
{
    int a,b;
    if(result!=0)  
    {
        a= (int) result;    printf("\na: %d\n",a);
        printf("\nFor a % b:\nEnter 'b' :\n");
        scanf("%d",&b);
    }

    else
    {
        printf("\n\nFor a % b:\nEnter 'a' & 'b' separated by space bar:\n");
        scanf("%d%d",&a,&b);
    }
    int remain = a%b;
    printf("(%d) % (%d) = ",a,b); return remain;
}
double power()
{
    double a,b;
    if(result!=0)  
    {
        a=result;   printf("\na: %lf\n",a);
        printf("\n\nFor a ^ b:\n\nEnter 'b' :\n");
        scanf("%lf",&b);
    }
    else
    {
        printf("\nFor a^b:\n\nEnter 'a' & 'b' separated by space bar (a>0,b>=0):\n");
        scanf("%lf%lf",&a,&b);
    }
    double exp = pow(a,b);
    if(a==0)
    {
        printf("\nInvalid Input for 'a' (a should be > 0)\n");
        return -1;
    }
    printf("\n(%lf) ^ (%lf) =",a,b); return exp;
}

int COLUMNS,rows;
void readmatrix(float arr[][COLUMNS])
{
    int col=COLUMNS;
    for(int i=1;i<=col-1;i++)
    {
        for(int j=1;j<=col;j++)
        {
			scanf("%f", &arr[i][j]);
        }
    }
    return;
}
void printmatrix(float arr[][COLUMNS])
{
    int col=COLUMNS;
    for(int i=1;i<=col-1;i++)
    {
        for(int j=1;j<=col;j++)
        {
            printf("%f\t",arr[i][j]);
        }
        printf("\n");
    }
    return;
}
void printsolution(float arr[][COLUMNS])
{
    for(int i=1;i<=rows;i++)
    {
        printf("\n x%d = %f ",i,(arr[i][COLUMNS]/arr[i][i]) );
    }
    return;
}
void checkforNoSol(float arr[][COLUMNS])
{
    int flag=0;
    for(int j=1;j<COLUMNS;j++)
    {
        for(int i=1;i<=COLUMNS;i++)
        {                                           
            if(i<j)                                 
            {                                       
                for(int k=j;k<COLUMNS;k++)          
                {                                   
                    if(arr[j][k]==0)                
                    flag=1;
                    else flag=0;
                }
            }
        }
    }
    if(flag==1)
    {
        printf("\nNo solution found for the Entered System of Linear Equations\n");
        printf("\n________________________________________________________________________\n");
        return;
    }

}
void checkforInfSol(float arr[][COLUMNS])
{
    int flag=0;
    for(int j=1;j<COLUMNS;j++)
    {
        for(int i=1;i<=COLUMNS;i++)
        {
            if(i<j)
            {
                for(int k=j;k<=COLUMNS;k++)
                {
                    if(arr[j][k]==0)
                    flag=1;
                    else flag=0;
                }
            }
        }
        if(flag==1) 
        {
            printf("\nInfinite solutions found for the Entered System of Linear Equations\n");
            printf("\n________________________________________________________________________\n");
            return;
        }
    }
    if(flag==1)
    {
        printf("\nInfinite solutions found for the Entered System of Linear Equations\n");
        printf("\n________________________________________________________________________\n");
        return;
    }

}
void diagonal(float arr[][COLUMNS])
{
    int i,j,k,col=COLUMNS;
    float ratio;
    for(j=1;j<col;j++)
    {
        for(i=1;i<col;i++)
        {
            if(j!=i)
            {
                ratio=(arr[i][j])/arr[j][j];
                for(k=1;k<=col;k++)
                {
                    arr[i][k] = (arr[i][k]) - (ratio*arr[j][k]);
                }
            }
        }
    }
    return;
}
void eqMain()
{
    float arr[20][20],solution[20];
    printf("\nEnter number of Linear Equations\n\n");
    scanf("%d",&rows);
    COLUMNS=rows+1;
    printf("\nEnter coefficients of equations (in the form of Ax + By + Cz (and so on..) = D):\n");
    readmatrix(arr);
    printf("\n______________________________________\n\nEntered Coefficients: \n");
    printmatrix(arr);
    diagonal(arr);
    printf("\n______________________________________\n\nSolution of the entered equations:\n");
    checkforInfSol(arr);
    checkforNoSol(arr);
    printsolution(arr); 
    printf("\n________________________________________________________________________\n");
}

double normalized(int num,int den)
{
    int arr[2],smaller; double ans,numer,deno;
    arr[0]=num; arr[1]=den;


    if(num==0 && den==0)
    {printf("NOT DEFINED\n");   exit(0);}
    if(num==0 && den!=0)
    {
        printf("0");
        return 0;
    }
    if(num!=0 && den==0)
    {printf("NOT DEFINED\n");   exit(0);}


    if((num>0 && den>0) || num<0 && den<0)
    {
        int absnumer=abs(num), absdeno=abs(den);
        smaller=(absnumer<=absdeno)?absnumer:absdeno;int i=1;
        for(;i<=smaller;i++)
        {
            if(arr[0]%i==0 && arr[1]%i==0)
            {
                arr[0]=arr[0]/i;    arr[1]=arr[1]/i;
                i=1;
            }

        }
        numer=abs(arr[0]); deno=abs(arr[1]); ans = (numer/deno); arr[0]=numer;arr[1]=deno;
        printf("%d/%d",arr[0],arr[1]);
    }
    else 
    {
        int absnum=abs(num),absdeno=abs(den);
        smaller=(absnum<absdeno)?absnum:absdeno;int i=1;
        for(;i<=smaller;i++)
        {
            if(arr[0]%i==0 && arr[1]%i==0)
            {
                arr[0]=arr[0]/i;    arr[1]=arr[1]/i;
                i=1;
            }

        }
        numer=abs(arr[0]); deno=abs(arr[1]); ans = (numer/deno)*(-1); arr[0]=numer;arr[1]=deno;
        printf("-(%d/%d)",arr[0],arr[1]);
    }
    return ans;
}
double muldiv(int num1, int den1, int num2, int den2)
{
    int mulnum=num1*num2; double ans;
    int mulden=den1*den2;
    ans = normalized(mulnum,mulden);
    return ans;
}
double addsub(int num1, int den1, int num2, int den2)
{
    int addnum=0,addden=1;double ans;
    addnum = (num1*den2)+(num2*den1);
    addden = den1*den2;
    ans = normalized(addnum,addden);
    return ans;
}
double fractCalciMain()
{
    int numer1,deno1,numer2,deno2,mod; double ans;
    printf("Enter 1st number (Numerator and Denominator separated by Space Bar): \n");
    scanf("%d%d",&numer1,&deno1);
    printf("\nInputed Number 1: %d/%d\n",numer1,deno1);
    printf("\nEnter 2nd number (Numerator and Denominator separated by Space Bar): \n");
    scanf("%d%d",&numer2,&deno2);
    printf("\nInputed Number 2: %d/%d\n",numer2,deno2);
    if(deno1==0||deno2==0)
    { printf("\nINVALID INPUT FOUND!\nEither 1 or both inputed numbers is invalid\n"); exit(0);}
    printf("\n______________________________________________________\n");
    printf("Enter Operation:\n\n 1 for Addition\n 2 for Substraction\n 3 for Multiplication\n 4 for Division\n\n");
    scanf("%d",&mod);
    switch(mod)
    {
        case 1: printf("\n(%d/%d) + (%d/%d) = ",numer1,deno1,numer2,deno2);
                ans = addsub(numer1,deno1,numer2,deno2); break;
        
        case 2: printf("\n(%d/%d) - (%d/%d) = ",numer1,deno1,numer2,deno2);
                numer2=-numer2;
                ans = addsub(numer1,deno1,numer2,deno2);numer2=-numer2;break;

        case 3: printf("\n(%d/%d) x (%d/%d) = ",numer1,deno1,numer2,deno2);
                ans = muldiv(numer1,deno1,numer2,deno2); break;

        case 4: printf("\n(%d/%d) / (%d/%d) = ",numer1,deno1,numer2,deno2);
                ans = muldiv(numer1,deno1,deno2,numer2);break;
    }
    printf("\n______________________________________________________\n");
    printf(" \n\t 1 - Show in decimal\n\t 2 - Continue\n\n");
    int response; scanf("%d",&response);
    if(response == 1)   printf("\n  = %lf",(ans));
    return ans;
}

int cofactor();
int n;
void readdeter(int mat[n][n],int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            scanf("%d",&mat[i][j]);
        }
        //printf("\n");
    }
    return;
}
void printdeter(int mat[n][n],int n)
{
   for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            printf("  %d  ",mat[i][j]);
        }
        printf("\n");
    } 
    return;
}
int determinant(int N, int mat[][N])
{
    int i,sign=1,det=0,detOfCoF = 0;
    if( N == 1 )
    {
        det = mat[0][0];
    }
    else if( N == 2 )
    {
        det = mat[1][1]*mat[0][0] - mat[1][0]*mat[0][1];
    }
    else
    {
        det=0;
        for(i=0;i<N;i++)
        {
            detOfCoF = cofactor(N,mat,i);
            det+=(sign*mat[0][i]*detOfCoF);
            sign*=-1;
        }
    }
    return det;
}
int cofactor(int n,int mat[][n],int col)
{
    int N=n-1, CoF[N][N];
    int i,j,CoFrow=0,CoFcol=0,detOfCoF=0;
    for(i=1;i<n;i++)
    {
        CoFcol=0;
        for(j=0;j<n;j++)
        {
            if( j!=col )
            {
                CoF[CoFrow][CoFcol] = mat[i][j];
                CoFcol++;
            }
        }
        CoFrow++;
    }
    detOfCoF = determinant(N,CoF);
    return detOfCoF;
}
void deterMain()
{
    printf("\nEnter Order of Determinant:\n\n");
    scanf("%d",&n);
    printf("\n________________________________________________________________________________________\n");
    int mat[n][n];
    printf("\nEnter Determinant in row-wise order (Separated by Space):\n\n");
    readdeter(mat,n);
    printf("________________________________________________________________________________________\n");
    printf("\nEntered Matrix:\n\n");
    printdeter(mat,n);
    printf("\n________________________________________________________________________________________\n\n");
    printf("Determinant of the Entered Matrix: \n\n");
    int det = determinant(n,mat);
    printf("%d",det);
    printf("\n________________________________________________________________________________________\n\n");
}

double vect1[2],vect2[2],vector[2];const double pi=3.14159;
void cartToPolar()
{
    double polar[2];
    polar[0] = sqrt(vector[0]*vector[0] + vector[1]*vector[1]);
    polar[1] = atan(vector[1]/vector[0]);
    polar[1]*=(180/pi);
    printf("\n(%lf,%lf) into (r,theta) = (%lf, %lf)\n",vector[0],vector[1],polar[0],polar[1]);
}
void polarToCart()
{
    double cart[2];
    double r = vector[0];
    double theta = (vector[1]*pi)/180;
    cart[0]= r*cos(theta);
    cart[1]= r*sin(theta);
    printf("\n(%lf,%lf) into (x,y) = (%lf,%lf)\n",vector[0],vector[1],cart[0],cart[1]);
}
void input2()
{
    printf("\nENTER VECTOR 1:\n");
    scanf("%lf%lf",&vect1[0],&vect1[1]);
    printf("-------VECTOR 1------\n x: %lf\n y: %lf\n\n",vect1[0],vect1[1]);
    printf("ENTER VECTOR 2:\n");
    scanf("%lf%lf",&vect2[0],&vect2[1]);
    printf("-------VECTOR 2-------\n a: %lf\n b: %lf\n",vect2[0],vect2[1]);
}
void input1()
{
    printf("\nENTER VECTOR:\n");
    scanf("%lf%lf",&vector[0],&vector[1]);
    printf("\nVECTOR\n x: %lf\n y: %lf\n\n",vector[0],vector[1]);
}
double scalarProduct()
{
    double sp = (vect1[0]*vect2[0])+(vect1[1]*vect2[1]);
    return sp;
}
double vectorProduct()
{
    double vp = (vect1[0]*vect2[1])-(vect2[0]*vect1[1]);
    return vp;
}
void addvect()
{
    double add[2];
    add[0]=vect1[0]+vect2[0];
    add[1]=vect1[1]+vect2[1];
    printf("\nAddition of Vector 1 (%lf,%lf) and Vector 2 (%lf,%lf) \n  =   (%lf,%lf)\n\n",vect1[0],vect1[1],vect2[0],vect2[1],add[0],add[1]);
}
void subvect()
{
    double sub[2];
    sub[0]=vect1[0]-vect2[0];
    sub[1]=vect1[1]-vect2[1];
    printf("\nSubstraction of Vector 1 (%lf,%lf) and Vector 2 (%lf,%lf) \n  =   (%lf,%lf)\n\n",vect1[0],vect1[1],vect2[0],vect2[1],sub[0],sub[1]);
}
void multiplyvect()
{
    int constant;printf("\nEnter Constant: ");
    scanf("%d",&constant);
    printf("\nMultiplication of vector (%lf,%lf) with a constant %d = (%lf,%lf)\n",vector[0],vector[1],constant,(constant*vector[0]),(constant*vector[1]));
}
double angle()
{
    double dotpro = scalarProduct();
    double magOfa = sqrt((vect1[0]*vect1[0]) + (vect1[1]*vect1[1]));
    double magOfb = sqrt((vect2[0]*vect2[0]) + (vect2[1]*vect2[1]));
    double ang = acos( dotpro / ((magOfa)*(magOfb)) );
    ang*=(180/pi); 
    return ang;
}
void vectorMain()
{
    printf("Enter Operation:\n\n\t 1 - Convert Coordinates from Cartesian to Polar Form\n\t 2 - Convert Coordinates from Polar to Cartesian Form\n\t 3 - Addition of 2 Vectors\n\t 4 - Substraction of 2 Vectors\n\t 5 - Multiplication of a Vector with a Constant\n\t 6 - Scalar Product of 2 Vectors\n\t 7 - Vector Product of 2 Vectors\n\t 8 - Angle Between 2 Vectors\n\n");
    int mod; scanf("%d",&mod);
    switch(mod)
    {
        case 1: printf("\n______________________________________\n\nTO CONVERT CARTESIAN FORM TO POLAR FORM :\n\nEnter a vector in Cartesian Form:\n");
                input1();
                cartToPolar();break;

        case 2: printf("\n______________________________________\n\nTO CONVERT POLAR FORM TO CARTESIAN FORM :\n\nEnter a vector in Polar Form:\n");
                input1();
                polarToCart();break;

        case 3: printf("\n______________________________________\n\nADDITION OF 2 VECTORS :\n\nEnter a vectors in Cartesian Form:\n");
                input2();
                addvect();
                break;

        case 4: printf("\n______________________________________\n\nSUBSTRACTION OF 2 VECTORS :\n\nEnter a vectors in Cartesian Form:\n");
                input2();
                subvect();
                break;

        case 5: printf("\n______________________________________\n\nMULTIPLICATION OF A VECTOR WITH A CONSTANT :\n\nEnter a vector in Cartesian Form:\n");
                input1();
                multiplyvect();
                break;

        case 6: printf("\n______________________________________\n\nSCALAR/DOT PRODUCT OF 2 VECTORS :\n\nEnter a vectors in Cartesian Form:\n");
                input2();   double scapro = scalarProduct();
                printf("\nSCALER (DOT) PRODUCT OF VECTOR 1 AND VECTOR 2: %lf\n",scapro);
                break;
        case 7: printf("\n______________________________________\n\nVECTOR/CROSS PRODUCT OF 2 VECTORS :\n\nEnter a vectors in Cartesian Form:\n");
                input2();   double vectpro=vectorProduct();
                printf("\nVECTOR (CROSS) PRODUCT OF VECTOR 1 AND VECTOR 2: %lf\n",vectpro);
                break;
        
        case 8: printf("\n______________________________________\n\nANGLE BETWEEN 2 GIVEN VECTORS :\n\nEnter a vectors in Cartesian Form:\n");
                input2();   double ang = angle();
                printf("\nANGLE BETWEEN VECTOR 1 (%lf,%lf) & VECTOR 2 (%lf,%lf) = %lf",vect1[0],vect1[1],vect2[0],vect2[1],ang);
                break;
    }
    printf("\n______________________________________________________\n");
}

double sinInv()
{
    double val,res;
    printf("\nEnter input value:\t");
    scanf("%lf",&val);
    res = asin(val);
    res = res*180/pi;
    printf("SIN inv( %lf ) = ",val);
    return res;
}
double cosInv()
{
    double val,res;
    printf("\nEnter input value:\t");
    scanf("%lf",&val);
    res = acos(val);
    res = res*180/pi;
    printf("COS inv( %lf ) = ",val);
    return res; 
}
double tantrigo(int val1,float val2)
{
    int val3=val1;
    while(val3>=360)
    {
        val3 = val1-360;
    }
    if(val3==90||val3==270)
    {   printf("TAN ( %d ) = Not Defined",val1);    return 1000;   }
    else
    {   printf("TAN ( %d ) = ",val1); return tan(val2); }
}
double tanInv()
{
    double val,res;
    printf("\nEnter input value:\t");
    scanf("%lf",&val);
    res = atan(val);
    res = res*180/pi;
    printf("TAN inv( %lf ) = ",val);
    return res;
}
double cosectrigo(int val1,float val2)
{
    int val3=val1;
    while(val3>=360)
    {
        val3 = val1-360;
    }
    double ans = 1/(sin(val2));
    if(val3==0||val3==180)
    {   printf("COSEC ( %d ) = Not Defined",val1);  return 1000; }
    else
    {   printf("COSEC ( %d ) = ",val1); return ans; }
}
double cosecInv()
{
    double val,res;
    printf("\nEnter input value:\t");
    scanf("%lf",&val);
    if((int)val == 0)
    {
        printf("\nCOSEC inv( %lf ) = Not Defined",val);
        return 1000;
    }
    res = asin(1/val);
    res = res*180/pi;
    printf("\nCOSEC inv( %lf ) = ",val);
    return res;
}
double sectrigo(int val1,float val2)
{
    int val3=val1;
    while(val3>=360)
    {
        val3 = val1-360;
    }
    double ans = 1/(cos(val2));
    if(val3==90||val3==270)
    {   printf("SEC ( %d ) = Not Defined",val1);    return 1000; }
    else
    {   printf("SEC ( %d ) = ",val1); return ans;   }
}
double secInv()
{
    double val,res;
    printf("\nEnter input value:\t");
    scanf("%lf",&val);
    if((int)val == 0)
    {
        printf("\nSEC inv( %lf ) = Not Defined",val);
        return 1000;
    }
    res = acos(1/val);
    res = res*180/pi;
    printf("\nSEC inv( %lf ) = ",val);
    return res;
}
double cottrigo(int val1,float val2)
{
    int val3 = val1;
    while(val3>=360)
    {
        val3 = val1-360;
    }
    double ans = (cos(val2)/sin(val2));
    if(val3==0||val3==180)
    {   printf("COT ( %d ) = Not Defined",val1);    return 1000;    }
    else
    {   printf("COT ( %d ) = ",val1); return ans;   }
}
double cotInv()
{
    double val,res;
    printf("\nEnter input value:\t");
    scanf("%lf",&val);
    if((int)val == 0)
    {
        printf("\nCOT inv( %lf ) = Not Defined",val);
        return 1000;
    }
    res = atan(1/val);
    res = res*180/pi;
    printf("\nCOT inv( %lf ) = ",val);
    return res;
}
double trigoMain()
{
    int mod,val1;float val2;double ans;
    printf("\nEnter Operation:\n\n\t 1 - SIN\t\t 2 - SIN (inv)\n\t 3 - COS\t\t 4 - COS (inv)\n\t 5 - TAN\t\t 6 - TAN (inv)\n\t 7 - COSEC\t\t 8 - COSEC (inv)\n\t 9 - SEC\t\t 10 - SEC (inv)\n\t 11 - COT\t\t 12 - COT (inv)\n\n");
    scanf("%d",&mod);
    if(mod==1||mod==3||mod==5||mod==7||mod==9||mod==11)
    {
        printf("\nEnter Input Value (in degrees)(in integers): \n");
        scanf("%d",&val1);
        val2=val1*(pi/180);
    }
    printf("\n______________________________________________________\n\n");
    switch(mod)
    {
        case 1: ans = sin(val2); printf("SIN ( %d ) = ",val1);   break;
        case 2: ans = sinInv();                                  break;
        
        case 3: ans = cos(val2); printf("COS ( %d ) = ",val1);   break;
        case 4: ans = cosInv(val1,val2);                         break;

        case 5: ans = tantrigo(val1,val2);                       break;
        case 6: ans = tanInv(val1,val2);                         break;  

        case 7: ans = cosectrigo(val1,val2);                     break;
        case 8: ans = cosecInv(val1,val2);                       break;
        
        case 9: ans = sectrigo(val1,val2);                       break;
        case 10: ans = secInv(val1,val2);                        break;

        case 11: ans = cottrigo(val1,val2);                      break;
        case 12: ans = cotInv(val1,val2);                        break;
    }
    return ans;
}

void programcontinue()
{
    printf("\nCurrent Memory: %lf\n",result);
    printf("\nEnter Further Action:\n\t C/c - All Clear (Clear Existing Memory) \n\t E/e - Exit\n\t Any other key to continue without Clearing Existing Memory\n");
    scanf("%s",&exitprog);
    if(exitprog =='C'|| exitprog =='c')   {   printf("\n______________________________________________________\n"); result = 0; counter = 0;   }
    else if(exitprog == 'E' || exitprog == 'e')
    {   
        printf("\n______________________________________________________________________________________________________________________\n");
        printf("\n\t\t\t\t\t  Calculator is Shutting Down\n\t\t\t\t\t\tG O O D B Y E ");
        printf("\n______________________________________________________________________________________________________________________\n");
        exit(0);
    }
    else {  printf("\n______________________________________________________\n"); return;  }
}

void main()
{
    int mode,op,modecheck;
    printf("\n______________________________________________________________________________________________________________________\n");
    printf("\n\t\t\t\t\t  Welcome to Priyansh's Scientific Calculator\n");
    printf("\n______________________________________________________________________________________________________________________\n");
    while(exitprog!='E'|| exitprog!='e')
    {
        counter++;
        printf("\nSelect Mode:\n\n\t 1 - Integer \n\t 2 - Fraction \n\t 3 - Linear Equation \n\t 4 - Determinant \n\t 5 - Complex / Vector \n\t 6 - Trigonometic\n\t 7 - Exit \n");
        if(counter>1){   printf("\n CAUTION: Changing the mode will make the Current Memory = 0\n\n"); }
        scanf("%d",&mode);
        if(counter > 1 && mode != modecheck)    result = 0;
        modecheck = mode;
        switch (mode)
        {
            case 1:
                printf("\n______________________________________________________\n");
                printf("Enter Operation:\n\n\t 1 - Addition\n\t 2 - Substraction\n\t 3 - Multiplication\n\t 4 - Division\n\t 5 - Remainder\n\t 6 - Power Calulations\n\t 7 - Back \n");
                scanf("%d",&op);
                printf("\n______________________________________________________\n");double result1 = result;
                switch(op)
                { 
                    case 1: result+= add();         printf("\n = %lf",result);          printf("\n______________________________________________________\n");   break;

                    case 2: result-= subs();        printf("\n = %lf",result);          printf("\n______________________________________________________\n");   break;

                    case 3: result*= multiply();    printf("\n = %lf",result);          printf("\n______________________________________________________\n");   break;

                    case 4: result/= divide();      printf("\n = %lf",result);          printf("\n______________________________________________________\n");   break;

                    case 5: result = mod();         printf(" %lf \n",result);           printf("\n______________________________________________________\n");   break;

                    case 6: result = power();
                            if(result!=-1)
                                printf(" %lf \n",result);
                                else
                                    result = 0;
                            printf("\n______________________________________________________\n");   break;

                    case 7: break;

                    default:printf("\nInvalid Operation Entry");                        printf("\n______________________________________________________\n");   break;
                }
                break;

            case 2:
                printf("\n______________________________________________________\n");
                result+= fractCalciMain();
                printf("\n______________________________________________________\n");
                break;

            case 3:
                printf("\n______________________________________________________\n");
                eqMain();
                break;

            case 4:
                printf("\n______________________________________________________\n");
                deterMain();
                break;

            case 5:
                printf("\n______________________________________________________\n");
                vectorMain();
                break;

            case 6:
                printf("\n______________________________________________________\n");
                double check = trigoMain(); int check2 = (int)check;    double resultTrigo;
                if(check2!=1000)
                {
                    resultTrigo = check;    result+=check;   printf(" %lf\n",resultTrigo);  printf("\n______________________________________________________\n");
                }
                else
                {
                    result+=0;
                    printf("\n______________________________________________________\n");
                }
                break;

            case 7:
                printf("\n______________________________________________________________________________________________________________________\n");
                printf("\n\t\t\t\t\t  Calculator is Shutting Down\n\t\t\t\t\t\tG O O D B Y E ");
                printf("\n______________________________________________________________________________________________________________________\n");
                exit(0);

            default:
                printf("\n______________________________________________________\n\nInvalid Mode Entry");
                printf("\n______________________________________________________\n");
                programcontinue();  exit(0);
        }
        programcontinue();
    }
    exit(0);
}