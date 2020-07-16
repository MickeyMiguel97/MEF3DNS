//All the local components needed to create the local Matrix M and local B

//Determinant - vil copia
float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

//A' - vil copia
void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

//B' - vil copia
void calculateBeta(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) =  1; 
    B.at(0).at(6) =  0; 
    B.at(0).at(7) =  0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) =  1; 
    B.at(0).at(10) =  0; 
    B.at(0).at(11) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0; 
    B.at(1).at(6) = 1;
    B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1;
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
    B.at(2).at(4) = -1;
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0;
    B.at(2).at(7) = 1;
    B.at(2).at(8) = -1;
    B.at(2).at(9) = 0;
    B.at(2).at(10) = 0;
    B.at(2).at(11) = 1;
}

//C' - vil copia
void calculateOmega(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0; C.at(0).at(3) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1; C.at(1).at(3) = 0;
    C.at(2).at(0) = -1; C.at(2).at(1) = 0; C.at(2).at(2) = 0; C.at(2).at(3) = 1;
}

//Jacobiano - vil copia
float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

//Special Matrix for my project
void calculateAries(int i,Matrix &Aries,mesh m){
    zeroes(Aries,12,3);
    element e = m.getElement(i);
    float cord1x = selectCoord(EQUIS,selectNode(1,e,m));
    float cord2x = selectCoord(EQUIS,selectNode(2,e,m));
    float cord3x = selectCoord(EQUIS,selectNode(3,e,m));
    float cord4x = selectCoord(EQUIS,selectNode(4,e,m));
    float cord1y = selectCoord(YE,selectNode(1,e,m));
    float cord2y = selectCoord(YE,selectNode(2,e,m));
    float cord3y = selectCoord(YE,selectNode(3,e,m));
    float cord4y = selectCoord(YE,selectNode(4,e,m));

    Aries.at(0).at(0) = (2*cord1x) + cord2x + cord3x + cord4x - (2*cord1y) - cord2y - cord3y - cord4y;
    Aries.at(0).at(1) = 0;
    Aries.at(0).at(2) = 0;
    Aries.at(1).at(0) = cord1x + (2*cord2x) + cord3x + cord4x - cord1y - (2*cord2y) - cord3y - cord4y;
    Aries.at(1).at(1) = 0;
    Aries.at(1).at(2) = 0;
    Aries.at(2).at(0) = cord1x + cord2x + (2*cord3x) + cord4x - cord1y - cord2y - (2*cord3y) - cord4y;
    Aries.at(2).at(1) = 0;
    Aries.at(2).at(2) = 0;
    Aries.at(3).at(0) = cord1x + cord2x + cord3x + (2*cord4x) - cord1y - cord2y - cord3y - (2*cord4y);
    Aries.at(3).at(1) = 0;
    Aries.at(3).at(2) = 0;

    Aries.at(4).at(0) = 0;
    Aries.at(4).at(1) = (2*cord1x) + cord2x + cord3x + cord4x - (2*cord1y) - cord2y - cord3y - cord4y;
    Aries.at(4).at(2) = 0;
    Aries.at(5).at(0) = 0;
    Aries.at(5).at(1) = cord1x + (2*cord2x) + cord3x + cord4x - cord1y - (2*cord2y) - cord3y - cord4y;
    Aries.at(5).at(2) = 0;
    Aries.at(6).at(0) = 0;
    Aries.at(6).at(1) = cord1x + cord2x + (2*cord3x) + cord4x - cord1y - cord2y - (2*cord3y) - cord4y;
    Aries.at(6).at(2) = 0;
    Aries.at(7).at(0) = 0;
    Aries.at(7).at(1) = cord1x + cord2x + cord3x + (2*cord4x) - cord1y - cord2y - cord3y - (2*cord4y);
    Aries.at(7).at(2) = 0;

    Aries.at(8).at(0) = 0;
    Aries.at(8).at(1) = 0;
    Aries.at(8).at(2) = (2*cord1x) + cord2x + cord3x + cord4x - (2*cord1y) - cord2y - cord3y - cord4y;
    Aries.at(9).at(0) = 0;
    Aries.at(9).at(1) = 0;
    Aries.at(9).at(2) = cord1x + (2*cord2x) + cord3x + cord4x - cord1y - (2*cord2y) - cord3y - cord4y;
    Aries.at(10).at(0) = 0;
    Aries.at(10).at(1) = 0;
    Aries.at(10).at(2) = cord1x + cord2x + (2*cord3x) + cord4x - cord1y - cord2y - (2*cord3y) - cord4y;
    Aries.at(11).at(0) = 0;
    Aries.at(11).at(1) = 0;
    Aries.at(11).at(2) = cord1x + cord2x + cord3x + (2*cord4x) - cord1y - cord2y - cord3y - (2*cord4y);
}

void calculateAcuario(int i,Matrix &Acuario,mesh m){
    zeroes(Acuario,12,3);
    element e = m.getElement(i);
    float cord1x = selectCoord(EQUIS,selectNode(1,e,m));
    float cord2x = selectCoord(EQUIS,selectNode(2,e,m));
    float cord3x = selectCoord(EQUIS,selectNode(3,e,m));
    float cord4x = selectCoord(EQUIS,selectNode(4,e,m));
    float cord1y = selectCoord(YE,selectNode(1,e,m));
    float cord2y = selectCoord(YE,selectNode(2,e,m));
    float cord3y = selectCoord(YE,selectNode(3,e,m));
    float cord4y = selectCoord(YE,selectNode(4,e,m));

    Acuario.at(0).at(0) = (2*cord1x) + cord2x + cord3x + cord4x + (2*cord1y) + cord2y + cord3y + cord4y;
    Acuario.at(0).at(1) = 0;
    Acuario.at(0).at(2) = 0;
    Acuario.at(1).at(0) = cord1x + (2*cord2x) + cord3x + cord4x + cord1y + (2*cord2y) + cord3y + cord4y;
    Acuario.at(1).at(1) = 0;
    Acuario.at(1).at(2) = 0;
    Acuario.at(2).at(0) = cord1x + cord2x + (2*cord3x) + cord4x + cord1y + cord2y + (2*cord3y) + cord4y;
    Acuario.at(2).at(1) = 0;
    Acuario.at(2).at(2) = 0;
    Acuario.at(3).at(0) = cord1x + cord2x + cord3x + (2*cord4x) + cord1y + cord2y + cord3y + (2*cord4y);
    Acuario.at(3).at(1) = 0;
    Acuario.at(3).at(2) = 0;

    Acuario.at(4).at(0) = 0;
    Acuario.at(4).at(1) = (2*cord1x) + cord2x + cord3x + cord4x + (2*cord1y) + cord2y + cord3y + cord4y;
    Acuario.at(4).at(2) = 0;
    Acuario.at(5).at(0) = 0;
    Acuario.at(5).at(1) = cord1x + (2*cord2x) + cord3x + cord4x + cord1y + (2*cord2y) + cord3y + cord4y;
    Acuario.at(5).at(2) = 0;
    Acuario.at(6).at(0) = 0;
    Acuario.at(6).at(1) = cord1x + cord2x + (2*cord3x) + cord4x + cord1y + cord2y + (2*cord3y) + cord4y;
    Acuario.at(6).at(2) = 0;
    Acuario.at(7).at(0) = 0;
    Acuario.at(7).at(1) = cord1x + cord2x + cord3x + (2*cord4x) + cord1y + cord2y + cord3y + (2*cord4y);
    Acuario.at(7).at(2) = 0;

    Acuario.at(8).at(0) = 0;
    Acuario.at(8).at(1) = 0;
    Acuario.at(8).at(2) = (2*cord1x) + cord2x + cord3x + cord4x + (2*cord1y) + cord2y + cord3y + cord4y;
    Acuario.at(9).at(0) = 0;
    Acuario.at(9).at(1) = 0;
    Acuario.at(9).at(2) = cord1x + (2*cord2x) + cord3x + cord4x + cord1y + (2*cord2y) + cord3y + cord4y;
    Acuario.at(10).at(0) = 0;
    Acuario.at(10).at(1) = 0;
    Acuario.at(10).at(2) = cord1x + cord2x + (2*cord3x) + cord4x + cord1y + cord2y + (2*cord3y) + cord4y;
    Acuario.at(11).at(0) = 0;
    Acuario.at(11).at(1) = 0;
    Acuario.at(11).at(2) = cord1x + cord2x + cord3x + (2*cord4x) + cord1y + cord2y + cord3y + (2*cord4y);
}

void calculatePiscis(int i,Matrix &Piscis,mesh m){
    zeroes(Piscis,3,12);
    element e = m.getElement(i);
    float cord1x = selectCoord(EQUIS,selectNode(1,e,m));
    float cord2x = selectCoord(EQUIS,selectNode(2,e,m));
    float cord3x = selectCoord(EQUIS,selectNode(3,e,m));
    float cord4x = selectCoord(EQUIS,selectNode(4,e,m));
    float cord1y = selectCoord(YE,selectNode(1,e,m));
    float cord2y = selectCoord(YE,selectNode(2,e,m));
    float cord3y = selectCoord(YE,selectNode(3,e,m));
    float cord4y = selectCoord(YE,selectNode(4,e,m));

    Piscis.at(0).at(0) = 3*pow(cord1x,2) + 2*cord1x*(cord2x+cord3x+cord4x) + pow(cord2x,2) + cord2x*(cord3x+cord4x) + pow(cord3x,2) + cord3x*cord4x + pow(cord4x,2) + 3*pow(cord1y,2) + 2*cord1y*(cord2y+cord3y+cord4y) + pow(cord2y,2) + cord2y*(cord3y+cord4y) + pow(cord3y,2) + cord3y*cord4y + pow(cord4y,2);
    Piscis.at(0).at(1) = pow(cord1x,2) + cord1x*(2*cord2x+cord3x+cord4x) + 3*pow(cord2x,2) + 2*cord2x*(cord3x+cord4x) + pow(cord3x,2) + cord3x*cord4x + pow(cord4x,2) + pow(cord1y,2) + cord1y*(2*cord2y+cord3y+cord4y) + 3*pow(cord2y,2) + 2*cord2y*(cord3y+cord4y) + pow(cord3y,2) + cord3y*cord4y + pow(cord4y,2);
    Piscis.at(0).at(2) = pow(cord1x,2) + cord1x*(cord2x+2*cord3x+cord4x) + pow(cord2x,2) + cord2x*(2*cord3x+cord4x) + 3*pow(cord3x,2) + 2*cord3x*cord4x + pow(cord4x,2) + pow(cord1y,2) + cord1y*(cord2y+2*cord3y+cord4y) + pow(cord2y,2) + cord2y*(2*cord3y+cord4y) + 3*pow(cord3y,2) + 2*cord3y*cord4y + pow(cord4y,2);
    Piscis.at(0).at(3) = pow(cord1x,2) + cord1x*(cord2x+cord3x+2*cord4x) + pow(cord2x,2) + cord2x*(cord3x+2*cord4x) + pow(cord3x,2) + cord3x*2*cord4x + 3*pow(cord4x,2) + pow(cord1y,2) + cord1y*(cord2y+cord3y+2*cord4y) + pow(cord2y,2) + cord2y*(cord3y+2*cord4y) + pow(cord3y,2) + cord3y*2*cord4y + 3*pow(cord4y,2);
    Piscis.at(0).at(4) = 0;
    Piscis.at(0).at(5) = 0;
    Piscis.at(0).at(6) = 0;
    Piscis.at(0).at(7) = 0;
    Piscis.at(0).at(8) = 0;
    Piscis.at(0).at(9) = 0;
    Piscis.at(0).at(10) = 0;
    Piscis.at(0).at(11) = 0;

    Piscis.at(1).at(0) = 0;
    Piscis.at(1).at(1) = 0;
    Piscis.at(1).at(2) = 0;
    Piscis.at(1).at(3) = 0;
    Piscis.at(1).at(4) = 3*pow(cord1x,2) + 2*cord1x*(cord2x+cord3x+cord4x) + pow(cord2x,2) + cord2x*(cord3x+cord4x) + pow(cord3x,2) + cord3x*cord4x + pow(cord4x,2) + 3*pow(cord1y,2) + 2*cord1y*(cord2y+cord3y+cord4y) + pow(cord2y,2) + cord2y*(cord3y+cord4y) + pow(cord3y,2) + cord3y*cord4y + pow(cord4y,2);
    Piscis.at(1).at(5) = pow(cord1x,2) + cord1x*(2*cord2x+cord3x+cord4x) + 3*pow(cord2x,2) + 2*cord2x*(cord3x+cord4x) + pow(cord3x,2) + cord3x*cord4x + pow(cord4x,2) + pow(cord1y,2) + cord1y*(2*cord2y+cord3y+cord4y) + 3*pow(cord2y,2) + 2*cord2y*(cord3y+cord4y) + pow(cord3y,2) + cord3y*cord4y + pow(cord4y,2);
    Piscis.at(1).at(6) = pow(cord1x,2) + cord1x*(cord2x+2*cord3x+cord4x) + pow(cord2x,2) + cord2x*(2*cord3x+cord4x) + 3*pow(cord3x,2) + 2*cord3x*cord4x + pow(cord4x,2) + pow(cord1y,2) + cord1y*(cord2y+2*cord3y+cord4y) + pow(cord2y,2) + cord2y*(2*cord3y+cord4y) + 3*pow(cord3y,2) + 2*cord3y*cord4y + pow(cord4y,2);
    Piscis.at(1).at(7) = pow(cord1x,2) + cord1x*(cord2x+cord3x+2*cord4x) + pow(cord2x,2) + cord2x*(cord3x+2*cord4x) + pow(cord3x,2) + cord3x*2*cord4x + 3*pow(cord4x,2) + pow(cord1y,2) + cord1y*(cord2y+cord3y+2*cord4y) + pow(cord2y,2) + cord2y*(cord3y+2*cord4y) + pow(cord3y,2) + cord3y*2*cord4y + 3*pow(cord4y,2);
    Piscis.at(1).at(8) = 0;
    Piscis.at(1).at(9) = 0;
    Piscis.at(1).at(10) = 0;
    Piscis.at(1).at(11) = 0;

    Piscis.at(2).at(0) = 0;
    Piscis.at(2).at(1) = 0;
    Piscis.at(2).at(2) = 0;
    Piscis.at(2).at(3) = 0;
    Piscis.at(2).at(4) = 0;
    Piscis.at(2).at(5) = 0;
    Piscis.at(2).at(6) = 0;
    Piscis.at(2).at(7) = 0;
    Piscis.at(2).at(8) = 3*pow(cord1x,2) + 2*cord1x*(cord2x+cord3x+cord4x) + pow(cord2x,2) + cord2x*(cord3x+cord4x) + pow(cord3x,2) + cord3x*cord4x + pow(cord4x,2) + 3*pow(cord1y,2) + 2*cord1y*(cord2y+cord3y+cord4y) + pow(cord2y,2) + cord2y*(cord3y+cord4y) + pow(cord3y,2) + cord3y*cord4y + pow(cord4y,2);
    Piscis.at(2).at(9) = pow(cord1x,2) + cord1x*(2*cord2x+cord3x+cord4x) + 3*pow(cord2x,2) + 2*cord2x*(cord3x+cord4x) + pow(cord3x,2) + cord3x*cord4x + pow(cord4x,2) + pow(cord1y,2) + cord1y*(2*cord2y+cord3y+cord4y) + 3*pow(cord2y,2) + 2*cord2y*(cord3y+cord4y) + pow(cord3y,2) + cord3y*cord4y + pow(cord4y,2);
    Piscis.at(2).at(10) = pow(cord1x,2) + cord1x*(cord2x+2*cord3x+cord4x) + pow(cord2x,2) + cord2x*(2*cord3x+cord4x) + 3*pow(cord3x,2) + 2*cord3x*cord4x + pow(cord4x,2) + pow(cord1y,2) + cord1y*(cord2y+2*cord3y+cord4y) + pow(cord2y,2) + cord2y*(2*cord3y+cord4y) + 3*pow(cord3y,2) + 2*cord3y*cord4y + pow(cord4y,2);
    Piscis.at(2).at(11) = pow(cord1x,2) + cord1x*(cord2x+cord3x+2*cord4x) + pow(cord2x,2) + cord2x*(cord3x+2*cord4x) + pow(cord3x,2) + cord3x*2*cord4x + 3*pow(cord4x,2) + pow(cord1y,2) + cord1y*(cord2y+cord3y+2*cord4y) + pow(cord2y,2) + cord2y*(cord3y+2*cord4y) + pow(cord3y,2) + cord3y*2*cord4y + 3*pow(cord4y,2);

}

void calculateLeo(int i,Matrix &Leo,mesh m){
    zeroes(Leo,12,3);

    Leo.at(0).at(0) = 1;
    Leo.at(0).at(1) = 1;
    Leo.at(0).at(2) = 1;
    Leo.at(1).at(0) = 1;
    Leo.at(1).at(1) = 1;
    Leo.at(1).at(2) = 1;
    Leo.at(2).at(0) = 1;
    Leo.at(2).at(1) = 1;
    Leo.at(2).at(2) = 1;
    Leo.at(3).at(0) = 1;
    Leo.at(3).at(1) = 1;
    Leo.at(3).at(2) = 1;

    Leo.at(4).at(0) = 1;
    Leo.at(4).at(1) = 1;
    Leo.at(4).at(2) = 1;
    Leo.at(5).at(0) = 1;
    Leo.at(5).at(1) = 1;
    Leo.at(5).at(2) = 1;
    Leo.at(6).at(0) = 1;
    Leo.at(6).at(1) = 1;
    Leo.at(6).at(2) = 1;
    Leo.at(7).at(0) = 1;
    Leo.at(7).at(1) = 1;
    Leo.at(7).at(2) = 1;

    Leo.at(8).at(0) = 1;
    Leo.at(8).at(1) = 1;
    Leo.at(8).at(2) = 1;
    Leo.at(9).at(0) = 1;
    Leo.at(9).at(1) = 1;
    Leo.at(9).at(2) = 1;
    Leo.at(10).at(0) = 1;
    Leo.at(10).at(1) = 1;
    Leo.at(10).at(2) = 1;
    Leo.at(11).at(0) = 1;
    Leo.at(11).at(1) = 1;
    Leo.at(11).at(2) = 1;

}

void calculateGeminis(int i,Vector &Geminis,mesh m){
    zeroes(Geminis,4);

    Geminis.at(0) = 1;
    Geminis.at(1) = 1;
    Geminis.at(2) = 1;
    Geminis.at(3) = 1;
}

Matrix createLocalM(int e,mesh &m){
    Matrix matrixA,matrixI,matrixL,matrixG,matrixD;
    float J,Determinant;
    
    // [ A-I  L+G ]
    // [  D   0 ]
    

    //Matrix A
    Matrix AriesMatrix, Alpha, Beta;

    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);

    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    float real_a = (float) (J)/(120*Determinant);
    calculateAries(e,AriesMatrix,m);
    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);
    productRealMatrix(real_a, productMatrixMatrix(AriesMatrix,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixA);

    //Matrix I
    Matrix alpha_t, beta_t;
    transpose(Alpha, alpha_t);
    transpose(Beta,beta_t);
    element ele = m.getElement(e);
    float cord1y = selectCoord(YE,selectNode(1,ele,m));
    float cord2y = selectCoord(YE,selectNode(2,ele,m));
    float cord3y = selectCoord(YE,selectNode(3,ele,m));
    float cord4y = selectCoord(YE,selectNode(4,ele,m));
    float cord1z = selectCoord(ZETA,selectNode(1,ele,m));
    float cord2z = selectCoord(ZETA,selectNode(2,ele,m));
    float cord3z = selectCoord(ZETA,selectNode(3,ele,m));
    float cord4z = selectCoord(ZETA,selectNode(4,ele,m));
    float Libra = cord1y*(2*cord1z+cord2z+cord3z+cord4z) + cord2y*(cord1z+2*cord2z+cord3z+cord4z) + cord3y*(cord1z+cord2z+2*cord3z+cord4z) + cord4y*(cord1z+cord2z+cord3z+2*cord4z);
    float real_Libra = (float) (J*Libra)/(120*Determinant*Determinant);
    //El real negativo por como quedo al aplicar el MEF
    productRealMatrix(-real_Libra, productMatrixMatrix(beta_t,productMatrixMatrix(alpha_t,productMatrixMatrix(Alpha,Beta,3,3,12),3,3,12),12,3,12),matrixI);

    //Matrix L
    Matrix omega;
    calculateOmega(omega);
    float cord1X = selectCoord(EQUIS,selectNode(1,ele,m));
    float cord2X = selectCoord(EQUIS,selectNode(2,ele,m));
    float cord3X = selectCoord(EQUIS,selectNode(3,ele,m));
    float cord4X = selectCoord(EQUIS,selectNode(4,ele,m));
    float cord1Z = selectCoord(ZETA,selectNode(1,ele,m));
    float cord2Z = selectCoord(ZETA,selectNode(2,ele,m));
    float cord3Z = selectCoord(ZETA,selectNode(3,ele,m));
    float cord4Z = selectCoord(ZETA,selectNode(4,ele,m));
    float Capricornio = cord1X*(2*cord1Z+cord2Z+cord3Z+cord4Z) + cord2X*(cord1Z+2*cord2Z+cord3Z+cord4Z) + cord3X*(cord1Z+cord2Z+2*cord3Z+cord4Z) + cord4X*(cord1Z+cord2Z+cord3Z+2*cord4Z);
    float real_Capricornio = (float) (J*Capricornio)/(120*Determinant*Determinant);
    productRealMatrix(real_Capricornio, productMatrixMatrix(beta_t,productMatrixMatrix(alpha_t,productMatrixMatrix(Alpha,omega,3,3,4),3,3,4),12,3,4),matrixL);

    //Matrix G
    Matrix accuario;
    calculateAcuario(e,accuario,m);
    float real_a = (float) J/(120*Determinant);
    productRealMatrix(real_a,productMatrixMatrix(accuario,productMatrixMatrix(Alpha,omega,3,3,4),12,3,4),matrixG);

    //Matrix D
    Matrix omega_t,piscis;
    transpose(omega,omega_t);
    calculatePiscis(e,piscis,m);
    float real_p = (float) J/(360*Determinant);
    productRealMatrix(real_p,productMatrixMatrix(omega_t,productMatrixMatrix(alpha_t,piscis,3,3,12),4,3,12),matrixD);

    //Matrix M
    Matrix M;
    zeroes(M,16);
    ubicarSubMatriz(M,0,11,0,11, sumMatrix(matrixA,matrixI,12,12));
    ubicarSubMatriz(M,0,11,12,15,sumMatrix(matrixL,matrixG,12,4));
    ubicarSubMatriz(M,12,15,0,11,matrixD);

    return M;
}

//Gravedad - vil copia
void calculateF(Vector &f, mesh &m){
    zeroes(f,3);

    f.at(0) += m.getParameter(EXTERNAL_FORCE_X);
    f.at(1) += m.getParameter(EXTERNAL_FORCE_Y);
    f.at(2) += m.getParameter(EXTERNAL_FORCE_Z);

}

Vector createLocalb(int e,mesh &m){
    float J;
    Vector b,b_aux,f,geminis,b_total,geminis_aux;
    Matrix leo;

    calculateF(f, m);
    calculateGeminis(e,geminis,m);

    calculateLeo(e,leo,m);

    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    //vector F 12x1
    zeroes(b_aux,12);
    productMatrixVector(leo,f,b_aux);
    productRealVector((500*J)/24,b_aux,b);

    //Vector H 4x1
    zeroes(geminis_aux,4);
    productRealVector((100*J)/24,geminis,geminis_aux);

    zeroes(b_total,16);

    joinTwoVectors(b,geminis_aux,b_total);
    
    return b_total;
}
