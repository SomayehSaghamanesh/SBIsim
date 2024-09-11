#ifndef PARSECHEMICALFORMULA_H
#define PARSECHEMICALFORMULA_H

#include "Materials.h"

#include <QRegularExpression>
#include <QMap>
#include <QRegularExpressionMatch>

struct ElementData {
    QString name;
    double atomicMass;
};

// Function to initialize element data
QMap<QString, ElementData> initializeElementData() {
    QMap<QString, ElementData> elementData;
    // Populate element data (atomic numbers and masses)
    // Example entries:
    elementData["H"] = {"Hydrogen", 1.008};
    elementData["He"] = {"Helium", 4.0026};
    // Add more elements as needed
    return elementData;
}

void ParseChemicalFormula(const QString& Name, QVector<int>& Z, QVector<double>& R) {
    Z.clear();
    R.clear();

    QString name = Name.trimmed();
    if (name.isEmpty()) return;

    // Initialize element data if not done already
    static QMap<QString, ElementData> elementData = initializeElementData();

    // Step 1: Check if Name is an element symbol
    if (name.length() <= 2) {
        auto it = elementData.find(name);
        if (it != elementData.end()) {
            Z.append(distance(elementData.begin(), it) + 1); // Atomic number (index)
            R.append(1.0); // Weight ratio
            return;
        }
    }

    // Step 2: Check if Name is an element name or known compound
    if (elementData.contains(name)) {
        Z.append(distance(elementData.begin(), elementData.find(name)) + 1); // Atomic number (index)
        R.append(1.0); // Weight ratio
        return;
    }

    // Step 3: Check if Name is a compound formula
    QRegularExpression compoundRegex("([A-Z][a-z]*)(\\(\\d*\\.?\\d+\\))?");
    QRegularExpressionMatchIterator matchIterator = compoundRegex.globalMatch(name);

    QVector<int> compoundZ;
    QVector<double> compoundR;

    while (matchIterator.hasNext()) {
        QRegularExpressionMatch match = matchIterator.next();
        QString elementName = match.captured(1);
        double ratio = 1.0;

        if (match.lastCapturedIndex() == 2) {
            QString ratioStr = match.captured(2);
            ratioStr.remove(QRegularExpression("[\\(\\)]")); // Remove parentheses
            ratio = ratioStr.toDouble();
        }

        if (elementData.contains(elementName)) {
            compoundZ.append(distance(elementData.begin(), elementData.find(elementName)) + 1); // Atomic number (index)
            compoundR.append(ratio); // Weight ratio
        }
    }

    // Normalize compound ratios
    double totalRatio = accumulate(compoundR.begin(), compoundR.end(), 0.0);
    if (totalRatio != 0.0) {
        for (double& ratio : compoundR) {
            ratio /= totalRatio;
        }
    }

    // Step 4: Check if Name is a mixture of elements and/or compounds
    if (!compoundZ.isEmpty()) {
        Z.append(compoundZ);
        R.append(compoundR);
        return;
    }


    QString Name = "aBcD"; // Example string

    // Create the regular expression and the replacement pattern
    QRegularExpression re("(\\D)([A-Z])");
    QString replacement = "\\1\\12";

    // Perform the replacement twice
    QString s1 = Name.replace(re, replacement);
    s1 = s1.replace(re, replacement);

    qDebug() << "Original: " << Name;
    qDebug() << "Transformed: " << s1;
}




#endif // PARSECHEMICALFORMULA_H
