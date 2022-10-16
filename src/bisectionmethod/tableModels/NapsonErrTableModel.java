package bisectionmethod.tableModels;

import javax.swing.table.DefaultTableModel;

public class NapsonErrTableModel extends DefaultTableModel {

    public NapsonErrTableModel(Object[][] data, String[] header) {
        super(data, header);
    }

    boolean[] canEdit = new boolean[] {
            false, false, false, false, false, false, false
    };

    public boolean isCellEditable(int rowIndex, int columnIndex) {
        return canEdit[columnIndex];
    }
}
