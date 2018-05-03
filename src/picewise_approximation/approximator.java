
package picewise_approximation;

//import java.util.*;
//import java.io.*;

public class approximator {
    private int n;                  //Количество вершин исходной выборки
    private double h;               //Шаг изменения V при пересчетах
    private double eps;             //Признак близости Вершины параболы и Tx
    private double[] Tx, Ty;        //T
    private double[] Vx, Vy;        //V
    private double[][] ABC;         //Коэффициенты параболы
    private double[][] l;           //Коэффициенты прямых
    private double[] Ux, Uy;            //Вершины парабол
    private double[] Uh;
    private double minH = 1e-4;      //Критерий уменьшения h
    private double maxIter = 10000;  //max число итераций
    public approximator() {
        System.out.println("Привет, я создался");
    }

    approximator(double[] X, double[] Y, int h, int eps) {
        n = X.length;
        Tx = X;
        Ty = Y;
        this.h = h;
        this.eps = eps;
        l = new double[3][n - 1];


        // уравнение прямой Li: Ai*y+Bi*x+Ci = 0, где Ai = L[i][0];
        for (int i = 0; i < n - 1; i++) {
            l[0][i] = (Tx[i + 1] - Tx[i]);
            l[1][i] = (Ty[i] - Ty[i + 1]);
            l[2][i] = (Ty[i] * Tx[i + 1] - Ty[i + 1] * Tx[i]);
        }
    }

    //Возможно указать любой коэффициент (какую часть длины отрезка Ti,Ti+1 нужно отсупить, чтобы выбрать )
    public void initV(double xh) {
        // xh = 0.5;
        Vx = new double[n - 1];
        Vy = new double[n - 1];

        Vx[0] = Tx[0];
        Vy[0] = Ty[0];
        Vx[n - 2] = Tx[n - 1];
        Vy[n - 2] = Ty[n - 1];

        for (int i = 1; i < n - 2; i++) {
            Vx[i] = Tx[i] + (Tx[i + 1] - Tx[i]) * xh;
            Vy[i] = -(l[1][i] * Vx[i] + l[2][i]) / l[0][i];
        }
    }

    public void initV() {
        initV(0.5);
    }

    public void VtoString() {
        for (int i = 0; i < n - 1; i++) {
            System.out.printf("%f\t%f\n", Vx[i], Vy[i]);
        }
    }

    public void initParabolas() {
        ABC = new double[3][n - 2];
        Ux = new double[n - 2];
        Uy = new double[n - 2];
        Uh = new double[n - 2];
        double denom;

        double Xu, Yu, Hu, Hu2;
        for (int i = 0; i < n - 2; i++) {
            denom = (Tx[i] - Tx[i + 1]) * (Tx[i + 2] - Tx[i]) * (Vx[i] - Vx[i + 1]);
            ABC[0][i] = (Tx[i] * (Ty[i + 1] - Ty[i + 2]) + Tx[i + 1] * (Ty[i + 2] - Ty[i]) + Tx[i + 2] * (Ty[i] - Ty[i + 1])) / 2 / denom;
            ABC[1][i] = (Vx[i] * (Tx[i] - Tx[i + 1]) * (Ty[i + 2] - Ty[i + 1]) - Vx[i + 1] * (Tx[i + 2] - Tx[i + 1]) * (Ty[i] - Ty[i + 1])) / denom;
            ABC[2][i] = ((Ty[i] * Tx[i + 1] - Ty[i + 1] * Tx[i]) * (Tx[i + 2] - Tx[i + 1]) * (Vx[i] - Vx[i + 1]) + Vx[i] * Vx[i] * (Ty[i] - Ty[i + 1]) * (Tx[i + 2] - Tx[i + 1]) + Vx[i] * Vx[i] * Ty[i + 2] *
                    (Tx[i + 1] - Tx[i])) / denom;

            Xu = (2 * ABC[0][i] * Tx[i + 1] - ABC[1][i] + Math.sqrt(Math.pow(2 * ABC[0][i] * Tx[i + 1] - ABC[1][i], 2) + 12 * ABC[0][i] * (Ty[i + 1] - ABC[2][i]))) / 6 / ABC[0][i];
            Yu = ABC[0][i] * Xu * Xu + ABC[1][i] * Xu + ABC[2][i];
            Hu = Math.sqrt(Math.pow((Tx[i + 1] - Xu), 2) - Math.pow((Ty[i + 1] - Yu), 2));

            System.out.printf("%f %f f\r\n",Math.pow(2 * ABC[0][i] * Tx[i + 1] - ABC[1][i], 2), 12*ABC[0][i] * (Ty[i + 1] - ABC[2][i]));
            Xu = (2 * ABC[0][i] * Tx[i + 1] - ABC[1][i] - Math.sqrt(Math.pow(2 * ABC[0][i] * Tx[i + 1] - ABC[1][i], 2) + 12 * ABC[0][i] * (Ty[i + 1] - ABC[2][i]))) / 6 / ABC[0][i];
            Yu = ABC[0][i] * Xu * Xu + ABC[1][i] * Xu + ABC[2][i];
            Hu2 = Math.sqrt(Math.pow((Tx[i + 1] - Xu), 2) - Math.pow((Ty[i + 1] - Yu), 2));
           // System.out.printf("%f %f %f\r\n",Xu,Yu,Hu2);

            if (Hu < Hu2) {
                Uh[i] = Hu;
            } else {
                Uh[i] = Hu2;
            }
        }
    }

    public void reCalcParabols(int i) {

        double denom, Xu, Yu, Hu, Hu2;
        denom = (Tx[i] - Tx[i + 1]) * (Tx[i + 2] - Tx[i]) * (Vx[i] - Vx[i + 1]);
        ABC[0][i] = (Tx[i] * (Ty[i + 1] - Ty[i + 2]) + Tx[i + 1] * (Ty[i + 2] - Ty[i]) + Tx[i + 2] * (Ty[i] - Ty[i + 1])) / 2 / denom;
        ABC[1][i] = (Vx[i] * (Tx[i] - Tx[i + 1]) * (Ty[i + 2] - Ty[i + 1]) - Vx[i + 1] * (Tx[i + 2] - Tx[i + 1]) * (Ty[i] - Ty[i + 1])) / denom;
        ABC[2][i] = ((Ty[i] * Tx[i + 1] - Ty[i + 1] * Tx[i]) * (Tx[i + 2] - Tx[i + 1]) * (Vx[i] - Vx[i + 1]) + Vx[i] * Vx[i] * (Ty[i] - Ty[i + 1]) * (Tx[i + 2] - Tx[i + 1]) + Vx[i] * Vx[i] * Ty[i + 2] *
                (Tx[i + 1] - Tx[i])) / denom;

        Xu = (2 * ABC[0][i] * Tx[i + 1] - ABC[1][i] + Math.sqrt(Math.pow(2 * ABC[0][i] * Tx[i] - ABC[1][i], 2) + 12 * ABC[0][i] * (Ty[i + 1] - ABC[2][i]))) / 6 / ABC[0][i];
        Yu = ABC[0][i] * Xu * Xu + ABC[1][i] * Xu + ABC[2][i];
        Hu = Math.sqrt(Math.pow((Tx[i + 1] - Xu), 2) - Math.pow((Ty[i + 1] - Yu), 2));

        Xu = (2 * ABC[0][i] * Tx[i + 1] - ABC[1][i] - Math.sqrt(Math.pow(2 * ABC[0][i] * Tx[i + 1] - ABC[1][i], 2) + 12 * ABC[0][i] * (Ty[i + 1] - ABC[2][i]))) / 6 / ABC[0][i];
        Yu = ABC[0][i + 1] * Xu * Xu + ABC[1][i + 1] * Xu + ABC[2][i];
        Hu2 = Math.sqrt(Math.pow((Tx[i + 1] - Xu), 2) - Math.pow((Ty[i + 1] - Yu), 2));

        if (Hu < Hu2) {
            Uh[i] = Hu;
        } else {
            Uh[i] = Hu2;
        }
    }

    public double getMaxH(boolean value, int exclude) {
        double maxH = Uh[0];
        int index = 0;

        for (int i = 0; i < Uh.length; i++) {
            if ((Uh[i] > maxH) && (i != exclude)) {
                maxH = Uh[i];
                index = i;
            }
        }
        if (value) {
            return maxH;
        } else {
            return index;
        }
    }

    public double calcE() {
        //Функция вычисления Е. Пусть будет sum(h)
        double sum = 0;
        for (int i=0;i<n-2;i++)
        {
            sum+=Uh[i];

        }
        return sum;
    }

    public void process() {
        initParabolas();
        int moving = 0;
        int tosee = n;
        int inThisH = 0;
        double prevE;
        int counter=0;
        while (calcE() > eps) {
            counter++;
            if (counter==maxIter ){
                break;
            }
            // Взять отрезок с максимальной погрешностью
            tosee = (int) getMaxH(false, tosee);

            //выбираем какую вершину двигать
            if (tosee == 0) {
                moving = 1;
            }
            if (tosee == n - 2) {
                moving = -1;
            }

            if (moving == 0) {
                if (Uh[tosee - 1] < Uh[tosee + 1]) {
                    moving = -1;
                } else {
                    moving = 1;
                }
            }
            //Запоминаем текущее положение дел
            prevE = calcE();

            if (inThisH == 2) {
                h = h / 2;
                if (h < minH) {
                    break;
                }
            }

            //двигаем
            Vx[tosee + moving] = Vx[tosee + moving] - h * moving;
            Vy[tosee + moving] = -(l[1][tosee + moving] * Vx[tosee + moving] + l[2][tosee + moving]) / l[0][tosee + moving];

            //пересчитать две параболы
            reCalcParabols(tosee);
            reCalcParabols(tosee + moving);

            if (calcE() >= prevE) {
                //Возвращаем на место прошлую вершину и пересчитываем ту параболу
                Vx[tosee + moving] = Vx[tosee + moving] + h * moving;
                Vy[tosee + moving] = -(l[1][tosee + moving] * Vx[tosee + moving] + l[2][tosee + moving]) / l[0][tosee + moving];
                reCalcParabols(tosee + moving);

                //Попробуем передвинуть другую вершину
                moving = -1 * moving;
                Vx[tosee + moving] = Vx[tosee + moving] - h * moving;
                Vy[tosee + moving] = -(l[1][tosee + moving] * Vx[tosee + moving] + l[2][tosee + moving]) / l[0][tosee + moving];

                //пересчитать две параболы
                reCalcParabols(tosee);
                reCalcParabols(tosee + moving);

                if (calcE() >= prevE) {
                    //Если и с этой не повезло возвращаем все как было
                    Vx[tosee + moving] = Vx[tosee + moving] + h * moving;
                    Vy[tosee + moving] = -(l[1][tosee + moving] * Vx[tosee + moving] + l[2][tosee + moving]) / l[0][tosee + moving];
                    reCalcParabols(tosee);
                    reCalcParabols(tosee + moving);

                    inThisH++;
                } else {
                    inThisH = 0;
                }
            } else {
                inThisH = 0;
            }//else
        }//while
        System.out.println(calcE());
    }//process
    public double getFunc(double p_x) {
        if (p_x >= Vx[0] && p_x <= Vx[n - 2]) {
            for (int i = 0; i < n - 2; i++) {
                if (p_x >= Vx[i] && p_x <= Vx[i + 1]) {
                    return (ABC[0][i] * p_x * p_x + ABC[1][i] * p_x + ABC[2][i]);
                }
            }
        }
        return 0;   //что то нужно вернуть, даже если непопал ни в один диапозон
    }//getFunc
}//class
