package cm;
import java.io.*;
import java.util.Scanner;


import static java.lang.Math.*;
public class Main {
    //вычисление евклидовой нормы матрицы
    static double euclideanNorm(double[][] matrix, int n){
        double norm1 = 0;
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                norm1 += matrix[i][j] * matrix[i][j];
            }
        }
        return sqrt(norm1);
    }

    //вычисление матричной нормы, подчиненной векторной l1
    static double m_norm(double[][] matrix, int n){
        double max = Double.MIN_VALUE;
        for(int i = 0; i < n; i++) {
            double s1 = 0;
            for (int j = 0; j < n; j++) {
                s1 += abs(matrix[i][j]);
            }
            if(s1 > max)
                max = s1;
        }
        return max;
    }

    //перестановка столбцов матрциы
    static void swapColumns(double[][] matrix, int c1, int c2, int n){
       for(int i = 0; i < n; i++){
           double t = matrix[i][c1];
           matrix[i][c1] = matrix[i][c2];
           matrix[i][c2] = t;
       }
    }

    //перестановка строк матрицы
    static void swapRows(double[][] matrix, int r1, int r2, int n){
        for(int i = 0; i < n; i++){
            double[] t = matrix[r1];
            matrix[r1] = matrix[r2];
            matrix[r2] = t;
        }
    }

    //решение системы с верхней треугольной матрицей
    static void solveSystem(double[][] s, double[] d, double[] x, double[] b, int n){
        double z[] = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += z[j] * s[j][i];
            }
            z[i] = (b[i] - sum) / s[i][i];
        }
        double y[] = new double[n];
        for (int i = 0; i < n; i++) {
            y[i] = z[i] * d[i];
        }

        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < n; j++) {
                sum += x[j] * s[i][j];
            }
            x[i] = (y[i] - sum) / s[i][i];
        }
    }

    //вычиление чисел обусловленности гильбертвоых матриц
    static double hilbert(){
        for(int size = 2; size <= 7; size++){
            double matrix[][] = new double[size][size];
            for(int i = 0; i < size; i++){
                for(int j = 0; j < size; j++){
                    matrix[i][j] = 1.0f / (i + j + 1);
                }
            }
            double s[][] = new double[size][size];
            double d[] = new double[size];
            squareRootDecomposition(matrix, s, d,  null, 0.00001f, null, size);
            double invMatrix[][] = new double[size][size];
            getInverseMatrix(s, invMatrix, null, d, size);

            System.out.println("Размер матрицы Гильберта = " + size);
            System.out.printf("число обусловленности c Евклидовой нормой = %.2f \n",euclideanNorm(matrix, size)*euclideanNorm(invMatrix, size));
            System.out.printf("число обусловленности c нормой, подчиненной векторной l1, = %.2f \n", m_norm(matrix, size) * m_norm(invMatrix, size));
        }
        return 0;
    }

    static boolean squareRootDecomposition(double[][] A, double[][] s, double[] d, double[] b, double e, int[] row, int n){
        boolean isDecompositionPossible = true;
        for(int i = 0; i < n; i++){
            //вычисляем диагональный элемент
            double sum = 0;
            for(int j = 0; j < i; j++){
                sum += s[j][i] * s[j][i] * d[j];
            }
            d[i] = signum(A[i][i] - sum);
            s[i][i] = (double)sqrt(abs(A[i][i] - sum));

            //если на диагонали оказался ноль
            if (abs(s[i][i]) < e){
                System.out.println("На диагонали матрицы s получен 0.");
                int k;
                boolean isAnotherExist = false;
                for(k = i + 1; k < n; k++)
                    if(abs(A[k][k] - A[i][i]) > e){
                        isAnotherExist = true;
                        break;
                    }
                if(isAnotherExist){
                    swapColumns(A, i, k, n);
                    swapRows(A, i, k, n);
                    swapColumns(s, i, k, n);
                    double t = b[i];
                    b[i] = b[k];
                    b[k] = t;
                    int t2 = row[i];
                    row[i] = row[k];
                    row[k] = t2;
                    System.out.println("Произведена перестановка строк и столбцов " + i + " и " + k);
                    //рестарт итерации
                    i--; continue;
                }
                else {
                    isDecompositionPossible = false;
                    System.out.println("Продолжить вычисления не является возможным.");
                }
            }
            //дальше по столбцам строки
            for(int k = i + 1; k < n; k++){
                double sum2 = 0;
                for(int j = 0; j < i; j++){
                    sum2 += s[j][i] * s[j][k]*d[j];
                }
                s[i][k] = (A[i][k] - sum2)/ (s[i][i]* d[i]);
            }
        }
        return isDecompositionPossible;
    }

    //получение обратной матрицы
    static void getInverseMatrix(double[][] s, double[][] inverseMatrix, int[] row, double d[], int n){
        for(int i = 0; i < n; i++){
            double right[] = new double[n];
            right[i] = 1;
            double column[] = new double[n];
            solveSystem(s, d, column, right, n);
            for(int j = 0; j < n; j++)
                if(row != null)
                    inverseMatrix[row[j]][row[i]] = column[j];
                else
                    inverseMatrix[j][i] = column[j];
        }
    }



    public static void main(String[] args) throws IOException {
        hilbert();
        OutputStream os = System.out;
        PrintStream fos = new PrintStream(new FileOutputStream("out.txt"));
        System.setOut(fos);
        boolean isDecompositionPossible = true;
        Scanner fin = new Scanner(new File("in.txt"));
        int n = fin.nextInt();

        final double e = 0.00001f;
        double[][] A = new double[n][n];
        double[][] originalMatrix = new double[n][n];
        for (int i = 0; i < n; i++) originalMatrix[i] = A[i].clone();
        double b[] = new double[n];
        double originalB[] = b.clone();
        double x[] = new double[n];

        //запоминаем перестановки строк
        //i - номер строки в текущей матрице,
        //row[i] - номер строки в исходной матрице
        int row[] = new int[n];
        for(int i = 0; i < n; i++)
            row[i] = i;

        //чтение матрицы из файла
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = fin.nextDouble();
            }
            b[i] = fin.nextDouble();
        }

        double d[] = new double[n];
        double s[][] = new double[n][n];

        isDecompositionPossible = squareRootDecomposition(A, s, d, b, e, row, n);

        if(isDecompositionPossible) {
            System.out.println("Матрица S: ");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    System.out.printf("%6.3f ", s[i][j]);
                }
                System.out.println();
            }
            System.out.println("Матрица D: ");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i == j)
                        System.out.printf("%2.0f ", d[i]);
                    else
                        System.out.printf("%2.0f ", 0.0);
                }
                System.out.println();
            }

            /*
                Решение системы
             */
            solveSystem(s, d, x, b, n);

            for (int i = 0; i < n; i++) {
                System.out.printf("x%d = %.3f   ", row[i] + 1, x[i]);
            }
            System.out.println();

            /*
                Нахождение определителя
             */
            double det = 1;
            for(int i = 0; i < n; i++){
                det *= s[i][i] * (s[i][i] * d[i]);
            }
            System.out.printf("detA = %.2f \n", det);

            /*
                Нахождение обратной матрицы
             */
            double inverseMatrix[][] = new double[n][n];
            getInverseMatrix(s, inverseMatrix, row, d, n);
            System.out.println("Обратная Матрица: ");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    System.out.printf("%7.3f ", inverseMatrix[i][j]);
                }
                System.out.println();
            }
        }
    }
}
