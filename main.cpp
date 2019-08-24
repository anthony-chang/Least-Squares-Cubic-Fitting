#include <iostream>
#include <bits/stdc++.h>

using namespace std;

/* N = number of data points
 * arrx and arry are arrays for the x and y data points respectively
 * V is the Vandermonde matrix
 * VT is the transpose of the Vandermonde matrix
 * VTV = VT * V
 *
 * adj is the adjoint matrix
 * inv is the inverse of VTV
 * VTy = VT * arry
 * ans is the coefficient matrix
 *
 * minVal is the minimum data point entered
 */
int N;
double arrx[1000], arry[1000], V[1000][4], VT[4][1000], VTV[4][4];
double adj[4][4], inv[4][4], VTy[4], ans[4];
double minVal;

//find the cofactor of at A[p][q]
void cofactor(double A[4][4], double temp[4][4], int p, int q, int n) {
    int i = 0, j = 0;
    for(int r = 0; r < n; r++) //loop through all of matrix A
        for(int c = 0; c < n; c++) {
            if(r!=p && c != q) { //copy everything from A to temp if it
                                //is not in row p and column q
                temp[i][j++] = A[r][c];
                if(j == n-1) { //if it is on the last column
                    j = 0; //reset the column
                    ++i; //move to next row
                }
            }
        }
}

//recursively find the determinant of A
double determinant(double A[4][4], int n) {
    double ret = 0; //the return value
    if(n==1) //base case, for a 1x1 matrix
        return A[0][0];
    double temp[4][4]; //stores cofactors
    int sign = 1; //swaps between plus or minus
    for(int i = 0; i < 4; i++) { //loop over first row
        cofactor(A, temp, 0, i, n); //get the cofactor of A[0][i]
        ret += sign * A[0][i] * determinant(temp, n-1); //add to determinant
        sign *= -1; //terms are added with alternating signs
    }
    return ret;
}
//get the adjoint of A
void adjoint(double A[4][4]) {
    int sign = 1; //stores whether it should be positive or negative
    double temp [4][4]; //store the cofactors
    for(int i = 0; i < 4; i++) { //loop through the matrix
        for(int j = 0; j < 4; j++) {
            cofactor(A, temp, i, j, 4); //get the cofactor of A[i][j]
            sign = ((i+j)%2==0)? 1:-1; //positive if row+column is even, negative otherwise
            adj[j][i] = sign*determinant(temp, 3); //also flips the matrix, swaps row & column
        }
    }
}

int main() {
    minVal = 0x3F3F3F; //temporary value

    /****************USER INPUT****************/
    printf("Enter the number of data points: ");
    scanf("%d", &N);
    //x data points
    printf("Enter the values of the x data points: ");
    for(int i = 0; i < N; i++) {
        scanf("%lf", &arrx[i]);
        if(abs(arrx[i]) < minVal) //keeps track of the minimum data value
            minVal = arrx[i];
    }
    //y data points
    printf("Enter the values of the y data points: ");
    for(int i = 0; i < N; i++)
        scanf("%lf", &arry[i]);

    /****************CALCULATIONS****************/
    //creating the Vandermonde matrix
    for(int i = 0; i < N; i++) { //N rows
        for (int j = 0; j < 3; j++) { //4 columns
            V[i][j] = pow(arrx[i], 3-j); //3rd order
            VT[j][i] = V[i][j]; //create the transpose of V by flipping row and columns
        }
        V[i][3] = 1; //constant
        VT[3][i] = 1; //constant
    }

    //multiplying V^T by V
    for(int i = 0; i < 4; i++) //loop through rows of V^T
        for(int j = 0; j < 4; j++) //loop through columns of V
            for(int k = 0; k < N; k++) //loop through columns of V^T
                VTV[i][j] += VT[i][k] * V[k][j];


    /****************OUTPUT****************/
    printf("\n~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~\n");
    //display the Vandermonde Matrix
    printf("\nVandermonde Matrix, V:\n");
    for (int i = 0; i < N; i++) { //loop through V
        for (int j = 0; j < 4; j++)
            printf("%lf ", V[i][j]);
        cout<<endl;
    }

    printf("\nVandermonde Transpose Matrix, V^T:\n");
    //display the transpose of the Vandermonde Matrix
    for (int i = 0; i < 4; i++) { //loop through VT
        for (int j = 0; j < N; j++)
            printf("%lf ", VT[i][j]);
        cout<<endl;
    }

    printf("\n(V^T)(V) Matrix: \n");
    //display the transpose of the Vandermonde Matrix times the Vandermonde Matrix
    for(int i = 0; i < 4; i++) { //loop through VTV
        for(int j = 0; j < 4; j++)
            printf("%lf ", VTV[i][j]);
        cout<<endl;
    }

    /****************BONUS QUESTION****************/
    printf("\n~~~~~~~~~~~~~~~~~BONUS~~~~~~~~~~~~~~~~~\n");
    double det = determinant(VTV, 4); //find the determinant of VTV
    if(det <= 0.1*abs(minVal)) { //if < 10% of minimum value, it is de facto 0
                                // prevents floating point error/rounding error
        printf("Inverse of (V^T * V) does not exist: determinant = 0");
        return 0; //end
    }
    adjoint(VTV); //find the adjoint of VTV
    for(int i = 0; i < 4; i++) //loop through the inv
        for(int j = 0; j < 4; j++)
            inv[i][j] = adj[i][j]/det; //inverse = adjoint/determinant

    //multiply V^T by y
    for(int i = 0; i < 4; i++) {//loop through rows of VT
        for(int k = 0; k < N; k++) //loop through the columns of VT
            VTy[i] += VT[i][k] * arry[k];
    }
    //multiple inv by V^T*y
    for(int i = 0; i < 4; i++) { //loop through the rows of inv
        for(int k = 0; k < 4; k++) //loop through the columns of inv
            ans[i] += inv[i][k] * VTy[k]; //ans = (VT*V)^-1 * (VTy)
    }

    /****************OUTPUT FOR BONUS QUESTION****************/
    printf("\nInverse of (V^T)(V) Matrix:\n");
    //display (V^T * V)^-1
    for(int i = 0; i < 4; i++) { //loop through inv
        for(int j = 0; j < 4; j++)
            printf("%lf ", inv[i][j]);
        cout<<endl;
    }
    printf("\n(V^T)*y Matrix:\n");
    //display V^T * y
    for(int i = 0; i < 4; i++) { //loop through VTy
        printf("%lf\n", VTy[i]);
    }
    printf("\nCoefficient Matrix:\n");
    //display (VT*V)^-1 * (VTy)
    for(int i = 0; i < 4; i++) { //loop through ans
        printf("%lf\n", ans[i]);
    }

    //print final function
    printf("\nTherefore, the obtained cubic function is:\n");
    //remove + sign for negative coefficients
    printf("%lfx^3%s%lfx^2%s%lfx%s%lf", ans[0], (ans[1]>=0 ? " +":" "), ans[1], (ans[2]>=0 ? " +":" "), ans[2], (ans[3]>=0? " +":" "), ans[3]);

    return 0;
}