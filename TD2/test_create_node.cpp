
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "retrieval.hpp"

const int dT = 10;

using std::cout;
using std::endl;

bool test_create_node(int T) {
    std::cout << "Testing the function create_node() on " 
    << T << " tests ...\t\t\n";

    int nb_errors_1 = 0;
    int nb_errors_2 = 0;
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<int> Ndis(0, 1000);

    for (int i = 0; i < T; i++) {
        int idx = Ndis(generator);

        node* n = create_node(idx);
        // Test if we created the node correctly
        if ((n->idx != idx) || (n->left != NULL) || (n->right != NULL)) {
            nb_errors_1++;
        }
        delete n;
    }

    for (int i = 0; i < T; i++) {
        node* nl = create_node(101);
        node* nr = create_node(1369);
        int c = Ndis(generator);
        int m = Ndis(generator);
        int idx = Ndis(generator);
        node* n = create_node(c, m, idx, nl, nr);
        // Test if we created the node correctly
        if ((n->c != c) || (n->m != m) || (n->idx != idx) ||
            (n->left != nl) || (n->right != nr)) {
            nb_errors_2++;
        }
        delete n;
        delete nl;
        delete nr;
    }

    cout << "   #errors = " << nb_errors_1 << " in the first create_node" << endl;
    cout << "   #errors = " << nb_errors_2 << " in the second create_node" << endl;
    bool success = ((nb_errors_1 == 0) && (nb_errors_2 == 0));
    std::cout << (success ? "[OK]" : "[NOK]") << std::endl;
    return success;
}

void print_help_msg() {
    std::cout << std::endl
              << "USAGE: To test your implementation of the create_node functions run the program as follows:" << std::endl
              << std::endl
              << " ./test_create_node <tests>" << std::endl
              << std::endl
              << "The parameter is optional." << std::endl
              << "<tests>      - a positive integer- ," << std::endl
              << "\t    the number of tests to perform "<< std::endl
              << "\t    default value is " << dT << std::endl
              << std::endl
              << std::endl;
}

int main(int argc, char* argv[]) {
    // Read in the arguments
    print_help_msg();

    int arg = 1;
    int T = (argc > arg) ? std::stoi(argv[arg]) : dT;

        test_create_node(T);

    return 0;
}
