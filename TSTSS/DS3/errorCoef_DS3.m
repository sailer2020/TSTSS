function err = errorCoef_DS3(Z,C)

err = sum(sum( abs(Z-C) )) / (size(Z,1)*size(Z,2));