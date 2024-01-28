typedef struct {
    int index;
} Body; 


int main() {
    int i, j;
    Body item[5]; 
    float forces[5][5];

    for (i=0; i<5; i++) {
        for (j=0; j<5; j++) {
            forces[i][j] = gravity_force(item[i], item[j]);
            if (i == j) continue;
        }
    }
}