package cm;

import java.io.*;
import java.util.Random;
import java.util.Scanner;


import static java.lang.Math.*;

public class Main {
    interface Expression{
        int method(double[][] A, double[] x, double e, double[] b, int n);
    }
    //норма l1
    static double vectorNorm(double[] v, int n) {
        double res = 0;
        for(int i = 0; i < n; i++){
            res += abs(v[i]);
        }
        return res;
    }

    static int jacobi(double[][] A, double[] x, double e, double[] b, int n){
        double[] new_x = new double[n];
        double[] dx = new double[n];
        int amountOfIterations = 0;
        double normOfDX = e + 1;
        while(normOfDX >= e) {
            amountOfIterations++;
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for(int j = 0; j < n; j++){
                    if(i != j){
                        sum += A[i][j] * x[j];
                    }
                }
                new_x[i] = - sum / A[i][i] + b[i] / A[i][i];
            }
            //вычичсление разницы между старым и новым решением
            // и сохранение нового решения
            for (int i = 0; i < n; i++) {
                dx[i] = new_x[i] - x[i];
                x[i] = new_x[i];
            }
            normOfDX = vectorNorm(dx, n);
        }
        return amountOfIterations;
    }

    static int seidel(double[][] A, double[] x, double e, double[] b, int n){
        double[] new_x = new double[n];
        double[] dx = new double[n];
        int amountOfIterations = 0;
        double normOfDX = e + 1;
        while(normOfDX >= e) {
            amountOfIterations++;
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for(int j = 0; j < n; j++){
                    if(i > j){
                        sum += A[i][j] * new_x[j];
                    }
                    else if(i < j){
                        sum += A[i][j] * x[j];
                    }
                }
                new_x[i] = - sum / A[i][i] + b[i] / A[i][i];

            }
            //вычичсление разницы между старым и новым решением
            // и сохранение нового решения
            for (int i = 0; i < n; i++) {
                dx[i] = new_x[i] - x[i];
                x[i] = new_x[i];
            }
            normOfDX = vectorNorm(dx, n);
        }
        return amountOfIterations;
    }

    //вычисление невзяки
    static double[] discrepancy(double[][] A, double[] x, double[] b, int n){
        double[] discrepancy = new double[n];
        for (int i = 0; i < n; i++){
            discrepancy[i] = 0;
            for(int j = 0; j < n; j++){
                discrepancy[i] += A[i][j] * x[j];
            }
            discrepancy[i] -= b[i];
        }
        return discrepancy;
    }

    //Выполняется ли диагональное преобладание
    static boolean checkDiagonalDomination(double[][] A, int n){
        boolean res = true;
        for(int i = 0; i < n; i++){
            double sum = 0;
            for(int j = 0; j < n; j++){
                if(i != j)
                    sum += A[i][j];
            }
            if(A[i][i] <= sum)
                res = false;
        }
        return res;
    }

    static void solve(Expression method, double[][] A, double[] x, double e, double[] b, int n, double[] x0){
        if(!checkDiagonalDomination(A, n))
            System.out.println("Диагональное преобладание не выполняется");
        if (n >= 0) System.arraycopy(x0, 0, x, 0, n);
        int amountOfIterations = method.method(A, x, e, b, n);
        for (int i = 0; i < n; i++){
            System.out.printf("x%d = %5.3f ", i + 1,  x[i]);
        }
        System.out.println();
        double[] discrepancy = discrepancy(A, x, b, n);
        System.out.print("Норма невзяки: " + vectorNorm(discrepancy, n) + "\n");
        System.out.println("Количество итераций: " + amountOfIterations);
    }




    public static void main(String[] args) throws FileNotFoundException {
        double e = 0.000001;
        double[][] A = null;
        double[] b = null;
        double[] x = null;
        double[] x0 = null;
        int n = 1;

        int option;
        boolean isMatrixCreated = false;
        System.out.printf("Стандартная погрешность: %f\n", e);
        do {
            System.out.println("1. Считать матрицу из файла");
            System.out.println("2. Сгенирировать матрицу c диагональным преобладанием");
            System.out.println("3. Решение методом Якоби");
            System.out.println("4. Решение методом Зейделя");
            System.out.println("5. Задать точность");
            System.out.println("6. Выбрать начальное приближение");
            System.out.println("7. Выход");
            Scanner scanner = new Scanner(System.in);
            option = scanner.nextInt();

            switch (option) {
                case 1 -> {
                    Scanner fin = new Scanner(new File("in.txt"));
                    n = fin.nextInt();
                    A = new double[n][n];
                    b = new double[n];
                    x0 = new double[n];
                    x = new double[n];
                    //чтение матрицы из файла
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            A[i][j] = fin.nextDouble();
                        }
                        b[i] = fin.nextDouble();
                    }
                    isMatrixCreated = true;
                }

                case 2 -> {
                    System.out.print("Введите размерность матрицы: ");

                    n = scanner.nextInt();

                    A = new double[n][n];
                    b = new double[n];
                    x0 = new double[n];
                    x = new double[n];

                    //генерация матрицы
                    Random random = new Random();
                    double sum = 0;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            if (i != j) {
                                A[i][j] = random.nextDouble();
                                sum += abs(A[i][j]);

                            }
                            A[i][i] = sum + abs(random.nextDouble());

                        }
                        b[i] = random.nextDouble()*100000;
                    }
                    isMatrixCreated =  true;

                    /*PrintWriter pw = new PrintWriter("generatedMatrix.txt");
                    pw.println(n);
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            pw.print(A[i][j] + " ");
                        }
                        pw.println(b[i]);
                    }
                    pw.close();*/
                }


                //решение методом Якоби
                case 3 -> {
                    if (isMatrixCreated) {
                        solve(Main::jacobi, A, x, e, b, n, x0);
                    }

                }
                //Решение методом Зейделя
                case 4 -> {
                    if (isMatrixCreated) {
                        solve(Main::seidel, A, x, e, b, n, x0);
                    }
                }

                case 5 -> {
                    System.out.print("Введите точность: ");
                    e = scanner.nextDouble();
                }
                case 6 -> {
                    if (isMatrixCreated) {
                        int internalOption;
                        System.out.println("1. Заполнить нулями");
                        System.out.println("2. Заполнить столбцом свободных членов");
                        System.out.println("3. Заполнить b[i]/A[i][i]");
                        internalOption = scanner.nextInt();
                        switch (internalOption) {
                            case 1: {
                                for (int i = 0; i < n; i++) {
                                    x0[i] = 0;
                                }
                            }break;
                            case 2: {
                                System.arraycopy(b, 0, x0, 0, n);
                            }break;
                            case 3:
                                for (int i = 0; i < n; i++) {
                                    x0[i] = b[i] / A[i][i];
                                }
                        }
                    }
                }


            }
        }while (option != 7);

    }
}
