
package bisectionmethod;
import java.util.ArrayList;

public class BisectionMethod {

   // tabel array list
    private ArrayList<Object[]> tbl = new ArrayList<Object[]>();

    public Object[][] newtonNapson(Double Es, ArrayList<Double> x, int i) {

        // Buat array untuk menampung eror, tipenya Object karena isinya berbagai macam
        // tipe data
        Object[] err = new Object[7];

        // index
        err[0] = i;

        // x[i]
        err[1] = x.get(i);

        // f(x)
        err[2] = f.__(x.get(i));

        // f'(x)
        err[3] = f.d(x.get(i));

        // x[i+1] = x[i]-(f(x)/f'(x))
        x.add(i + 1, x.get(i) - ((Double) f.__(x.get(i)) / f.d(x.get(i))));

        // x[i+1]
        err[4] = x.get(i + 1);

        // ganti index
        i++;

        // hentikan iterasi jika masih pertama kali
        if (i == 1) {

            // simpan rekam eror ke dalam tabel
            tbl.add(err);

            // jalankan ulang
            return newtonNapson(Es, x, i);
        }

        // |Ea|=|(x(i)-x(i-1)) /x(i)| *100
        Double Ea = (Double) Math.abs((Double) (x.get(i) - x.get(i - 1)) / x.get(i)) * 100;

        // Ea
        err[5] = Ea;

        // simpan rekam error ke dalam tabel
        tbl.add(err);

        // lanjut iterasi hingga Ea < Es
        if (Ea > Es) {
            // status
            err[6] = "next iterration";

            // jalankan ulang
            return newtonNapson(Es, x, i);
            // karena ada return, maka baris berikutnya tidak ada berjalan
        }

        // jika Ea < Es, maka iterasi akan diselesaikan
        err[6] = "stop iterration";

        // ubah arraylist menjadi array biasa
        Object[][] tbl_arr = new Object[tbl.size()][err.length];
        for (int j = 0; j < tbl_arr.length; j++) {
            tbl_arr[j] = tbl.get(j);
        }

        return tbl_arr;
    }

}

/**
 * ini adalah class untuk menampung f(x) dan f'(x)
 */
class f {

    public static Double __(Double x) {
        return x * Math.exp(-1 * x) + 1;
    }

    public static Double d(Double x) {
        return (1 - x) * Math.exp(-1 * x);
    }
}
