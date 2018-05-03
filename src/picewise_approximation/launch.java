package picewise_approximation;

public class launch {
    public static void main (String [] args){
        int n=5;
        double[][] xy=new double [2][n];
       /* for(int i=0;i<n;i++){
            xy[0][i]=i+1;
            xy[1][i]=i+1;
        }*/
        xy[0][0]=0;xy[1][0]=0;
        xy[0][1]=2;xy[1][1]=1;
        xy[0][2]=3;xy[1][2]=-1;
        xy[0][3]=3.5;xy[1][3]=0;
        xy[0][4]=4;xy[1][4]=-2;

        approximator ap = new approximator(xy[0],xy[1],0,0);
        ap.initV();
        ap.process();
        System.out.println("Где то в формулах ошибка");
        //ap.VtoString();
        for (double x=xy[0][0]; x<=xy[0][n-1];x=x+0.2)
        {
   //         System.out.printf("(%f;%f) ",x,ap.getFunc(x));
        }
    }
}
