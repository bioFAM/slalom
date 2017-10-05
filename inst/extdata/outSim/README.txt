HDF5 files contain:
    data_file['Y'] = Y - observed data matrix (input)
    data_file['Pi'] = Pi - prior for weight being on (observed annotations) (input)
    data_file['X'] = X - factor (output)
    data_file['W'] = W - weight (output)
    if FNper>0:
        data_file['FN'] = FN - false negatives
        data_file['FP'] = FP
        data_file['TN'] = TN
        data_file['TP'] = TP
    data_file['Ion'] = Ion -  True annotations
    if Khidden>0:
        data_file['Xhidden'] = Xhidden - unannotated factors, weights and indicator which weights are active (output)
        data_file['WHidden'] = Whidden
        data_file['IonHidden'] = IonHidden
        data_file['Nhidden'] = Khidden - number of unannotated factors
        data_file['IonTerm'] = IonTerm - which factors are active (5 active, 15 inactive)
