package com.company;
import java.io.*;
import java.util.Locale;
import java.util.Scanner;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

public class Main {
    public static void main(String[] args) throws IOException {
        OutputStream os = System.out;
        PrintStream fos = new PrintStream(new FileOutputStream("out.txt"));
        System.setOut(fos);
        boolean isMatrixDegenerate = false;
        Scanner fin = new Scanner(new File("in.txt"));
        int n = fin.nextInt();

        final float e = 0.00001f;
        float[][] matrix = new float[n][n];
        float[][] originalMatrix = new float[n][n];
        for(int i = 0; i < n; i++) originalMatrix[i] = matrix[i].clone();
        float b[] = new float[n];
        float originalB[] = b.clone();
        float x[] = new float[n];
        int detSign = 1;
        fin.useLocale(Locale.US);
        //чтение матрицы из файла
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = fin.nextFloat();
            }
            b[i] = fin.nextFloat();
        }

        //запоминаем перестановки строк
        //i - номер строки в текущей матрице,
        //row[i] - номер строки в исходной матрице
        int row[] = new int[n];
        for(int i = 0; i < n; i++)
            row[i] = i;


        /*
            прямой ход
         */

        //приведение к треугольному виду
        //по строкам
        for (int i = 0; i < n - 1; i++) {
            //поиск главного элемента
            float max = matrix[i][i];
            int imax = i;
            for (int j = i; j < n; j++)
                if (abs(matrix[j][i]) > max) {
                    max = abs(matrix[j][i]);
                    imax = j;
                }
            //если строчка с максимальным элементом не текущая
            if (imax != i) {
                float[] t1 = matrix[imax];
                matrix[imax] = matrix[i];
                matrix[i] = t1;
                float t2 = b[imax];
                b[imax]= b[i];
                b[i] = t2;
                int t3 = row[i];
                row[i] = row[imax];
                row[imax] = t3;
                detSign *= -1;
            }
            //число на главной диагонали и числа ниже - нули
            //матрица вырождена, определитель равен 0
            if(abs(matrix[i][i]) < e){
                isMatrixDegenerate = true;
                break;
            }


            //по оставшимся строкам
            for (int j = i + 1; j < n; j++) {
                float koef = matrix[j][i] / matrix[i][i];
                //сохраняем значения коэффициента для вычисления обратной матрицы
                matrix[j][i] = koef;
                //по столбцам
                for (int k = i + 1; k < n; k++) {
                    matrix[j][k] -= matrix[i][k] * koef;
                }
                b[j] -= b[i] * koef;

            }
        }
        //правый нижний элемент оказался равен нулю
        if(matrix[n-1][n-1] == 0)
            isMatrixDegenerate = true;


        //вывод получившейся матрицы
        System.out.println("Приведенная матрица: ");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if(i <= j)
                    System.out.printf("%7.3f ", matrix[i][j]);
                else
                    System.out.printf("%7.3f ", 0.0);
            }
            System.out.printf("|%7.3f\n", b[i]);
        }


        if(!isMatrixDegenerate) {
            /*
                Обратный ход
             */
            for (int i = n - 1; i >= 0; --i) {

                float tempb = b[i];
                for (int j = n - 1; j > i; --j) {
                    tempb -= matrix[i][j] * x[j];
                }
                x[i] = tempb / matrix[i][i];
            }

            System.out.println("Решение:");
            for (float __ : x) {
                System.out.printf("%7.3f ", __);
            }
            System.out.println();
             /*
                Вычисление невязки
                Подставляем х в исходную матрицу
             */
            float discrepancy[] = new float[n];
            for(int i = 0; i < n; i++){
                float left = 0;
                for(int j = 0; j < n; j++){
                    left += originalMatrix[i][j] * x[j];
                }
                discrepancy[i] = originalB[i] - left;
            }

            System.out.println("Невязка: ");
            for(float __: discrepancy){
                System.out.printf("%7.3f ", __);
            }
            System.out.println();

            //вычисление определителя
            double det = 1;
            for (int i = 0; i < n; i++){
                det *= matrix[i][i];
            }
            det *= detSign;
            System.out.printf("detA = %.3f\n", det);

            /*
                Вычисление обратной матрицы
             */

            float inverseMatrix[][] = new float[n][n];
            for(int q = 0; q < n; q++){
                float newB[] = new float[n];
                newB[q] = 1;
                //исправляем правую часть
                //начинаем с q, т.к. выше в правой части нули
                for(int j = q; j < n - 1; j++){
                    for(int i = j + 1; i < n; i++){
                        newB[i] -= newB[j] * matrix[i][j];
                    }
                }

                /*
                    Обратный ход
                */
                float[] tempX = new float[n];
                for (int i = n - 1; i >= 0; --i) {
                    float tempb = newB[i];
                    for (int j = n - 1; j > i; --j) {
                        tempb -= matrix[i][j] * tempX[j];
                    }
                    tempX[i] = tempb / matrix[i][i];
                }
                for(int i = 0; i < n; i++)
                    inverseMatrix[i][row[q]] = tempX[i];

            }
            //вывод обратной матрицы
            System.out.println("Обратная матрица:");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    System.out.printf("%7.3f ", inverseMatrix[i][j]);
                }
                System.out.println();
            }

            //Вычисление евклидовых норм
            double norm1 = 0;
            double norm2 = 0;
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    norm1 += matrix[i][j] * matrix[i][j];
                    norm2 += inverseMatrix[i][j] * inverseMatrix[i][j];
                }
            }
            norm1 = sqrt(norm1);
            norm2 = sqrt(norm2);
            System.out.printf("||A|| * ||invA|| = %.3f", norm1 * norm2);
        }
        else {
            System.out.println("Матрица вырождена, detA = 0");
        }

    }
}
