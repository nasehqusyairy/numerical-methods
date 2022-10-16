public class GaussJordan {

  public boolean isGaussed = false;
  public boolean isJordaned = false;

  private double[][] matrixA = {
      { 3, -0.1, -0.2 },
      { 0.1, 7, -0.3 },
      { 0.3, -0.2, 10 }
  };

  private double[][] matrixB = {
      { 7.85 },
      { -19.3 },
      { 71.4 }
  };

  public int length() {
    return matrixB.length;
  }

  public boolean isSquare() {
    return (matrixA[0].length == matrixA.length);
  }

  public boolean isNonHomogen() {
    for (int i = 0; i < matrixA.length; i++) {
      if (matrixB[i][0] != 0) {
        return true;
      }
    }
    return false;
  }

  public void gauss() {
    /*
     * Berdasarkan pengamatan :
     * 
     * index diagonal => index baris = iterasi ke-i, index kolom = itersi-ke i
     * 
     * jumlah iterasi = panjang matrix - 1
     * 
     * baris yang terpengaruh pada setiap iterasi adalah baris yang berada di bawah
     * diagonal
     * 
     * jumlah k = jumah baris yang terpengaruh
     * 
     * jumlah k maksimum = jumlah iterasi
     * 
     */

    try {

      // cek persyaratan
      if (!isSquare()) {
        throw new Exception("Matrix tidak persegi");

      } else if (!isNonHomogen()) {
        throw new Exception("Matrix tidak non-Homogen");

      }
      // else if (isGaussed) {
      // throw new Exception("Matrix sudah pernah melalui operasi Gauss");

      // }

      // diketahui
      int panjang_matrix = matrixA.length;
      int jumlah_iterasi = panjang_matrix - 1;
      double[] k = new double[jumlah_iterasi];

      for (int i = 1; i <= jumlah_iterasi; i++) {

        // nilai i untuk index diagonal dikurangi 1 karena dimulai dari nol
        double diagonal = matrixA[i - 1][i - 1];

        int jumlah_baris_yang_akan_terpengaruh = panjang_matrix - i;

        // buat nilai k
        for (int m = 0; m < jumlah_baris_yang_akan_terpengaruh; m++) {

          // index baris = m+i karena eksekusi baris dimulai dari index ke-i
          k[m] = matrixA[m + i][i - 1] / diagonal;

        }

        // matrix a
        // eksekusi kolom
        for (int l = 0; l < panjang_matrix; l++) {

          // eksekusi baris
          for (int j = 0; j < jumlah_baris_yang_akan_terpengaruh; j++) {

            // index baris = j+i karena eksekusi baris dimulai dari index ke-i
            matrixA[j + i][l] -= k[j] * matrixA[i - 1][l];

          }

        }

        // matrix b
        // eksekusi baris
        for (int j = 0; j < jumlah_baris_yang_akan_terpengaruh; j++) {

          // index baris = j+i karena eksekusi baris dimulai dari index ke-i
          matrixB[j + i][0] -= k[j] * matrixB[i - 1][0];

        }

      }

      isGaussed = true;

      // jika tidak memenuhi persyaratan, jalankan catch
    } catch (Exception e) {
      System.err.println("\nError:\n" + e.getMessage() + "operasi Gauss tidak dapat dilakukan");
    }

  }

  public void jordan() {
    if (isGaussed) {
      // diketahui
      int panjang_matrix = matrixA.length;
      int jumlah_iterasi = panjang_matrix - 1;
      double[] k = new double[jumlah_iterasi];

      for (int i = jumlah_iterasi; i >= 1; i--) {

        int jumlah_baris_yang_akan_terpengaruh = i;

        // buat nilai k
        for (int m = 0; m < jumlah_baris_yang_akan_terpengaruh; m++) {

          k[m] = matrixA[m][i] / matrixA[m + 1][i];

        }

        // matrix a
        // eksekusi kolom
        for (int l = jumlah_iterasi; l >= 0; l--) {

          // eksekusi baris
          for (int j = 0; j < jumlah_baris_yang_akan_terpengaruh; j++) {

            matrixA[j][l] -= k[j] * matrixA[j + 1][l];

          }

        }

        // matrix b
        // eksekusi baris
        for (int j = 0; j < jumlah_baris_yang_akan_terpengaruh; j++) {

          matrixB[j][0] -= k[j] * matrixB[j + 1][0];

        }

      }
      sederhanakanJordan();
      isJordaned = true;
    } else {
      gauss();
      jordan();
    }
  }

  private void sederhanakanJordan() {
    for (int i = 0; i < matrixA.length; i++) {
      double diagonal = matrixA[i][i];

      for (int j = 0; j < matrixA.length; j++) {
        matrixA[i][j] *= 1 / diagonal;
      }

      matrixB[i][0] *= 1 / diagonal;

    }
  }

  // method untuk mencari kofaktor
  static void kofaktor(double mat[][], double temp[][], int p, int q, int n) {
    int i = 0, j = 0;

    for (int row = 0; row < n; row++) {
      for (int col = 0; col < n; col++) {

        // salin ke variabel sementara
        if (row != p && col != q) {
          temp[i][j++] = mat[row][col];

          // tambahkan index baris dan mulai ulang undex kolom
          if (j == n - 1) {
            j = 0;
            i++;
          }
        }
      }
    }
  }

  // method untuk mencari determinan matrix secara rekursif
  static double determinan(double mat[][], int n) {
    double D = 0;

    // apabila ukura matrix = 1
    if (n == 1)
      return mat[0][0];

    // variabel kofaktor
    double temp[][] = new double[n][n];

    // tanda positif/negatif
    int sign = 1;

    // Iterate untuk setiap elemen di baris pertama
    for (int f = 0; f < n; f++) {

      // menentukan kofaktor mat[0][f]
      kofaktor(mat, temp, 0, f, n);
      D += sign * mat[0][f] * determinan(temp, n - 1);

      // balikkan tanda
      sign = -sign;
    }

    return D;
  }

  // method untuk memulai determinan
  public double determinan() {
    return determinan(matrixA, matrixA.length);
  }

  public double subtitusiGauss(int i) {

    if (isJordaned) {
      return matrixB[i][0];
    }

    if (determinan() == 0) {
      System.err.println("\nError : Determinan matrix sama dengan nol");
      return 0;
    }

    double subs = matrixB[i][0];

    if (i == matrixB.length - 1) {
      return matrixB[i][0] / matrixA[i][i];
    }

    for (int j = i + 1; j < matrixA.length; j++) {
      subs -= matrixA[i][j] * subtitusiGauss(j);

    }

    return subs / matrixA[i][i];
  }

  // method untuk menampilkan hasil
  public void display() {
    System.out.println("MatrixA");
    System.out.println("===============");
    for (int i = 0; i < matrixA.length; i++) {
      for (double ds : matrixA[i]) {
        System.out.print(ds + " ");
      }
      System.out.println();
    }
    System.out.println();
    System.out.println("MatrixB");
    System.out.println("===============");
    for (int i = 0; i < matrixB.length; i++) {
      for (double ds : matrixB[i]) {
        System.out.print(ds + " ");
      }
      System.out.println();
    }
  }

  public static void main(String[] args) {
    GaussJordan obj = new GaussJordan();
    // System.out.println(obj.isNonHomogen());
    System.out.println("# Sebelum operasi:");
    obj.display();

    obj.gauss();
    System.out.println();

    System.out.println("# Sesudah operasi Gauss:");
    obj.display();
    System.out.println();

    System.out.println("# Determinan : " + obj.determinan());
    System.out.println();

    obj.jordan();
    System.out.println();

    System.out.println("# Sesudah operasi Gauss-Jordan:");
    obj.display();
    System.out.println();

    for (int i = 0; i < obj.length(); i++) {
      System.out.println("x" + (i + 1) + " = " + obj.subtitusiGauss(i));
    }

  }

}
