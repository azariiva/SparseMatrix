#include "RBTree.hpp"
#include <cstddef>
#include <iostream>

int main() {
    RBTree<size_t, int> tree = RBTree<size_t, int>();
    tree.insert(10, 10);
    tree.insert(3, 3);
    tree.insert(7, 7);
    tree.insert(4, 4);
    tree.insert(20, 20);
    tree.insert(15, 15);
    tree.insert(11, 11);
    tree.insert(40, 40);
    tree.run();
    std::cout << std::endl;

    tree.remove(11);
    tree.remove(3);
    tree.remove(4);
    tree.run();
    std::cout << std::endl;

    tree.remove(15);
    tree.remove(7);
    tree.remove(20);
    tree.run();
    std::cout << std::endl;

    tree = RBTree<size_t, int>();
    tree.insert(10, 10);
    tree.insert(11, 11);
    tree.insert(12, 12);
    tree.insert(13, 13);
    tree.insert(14, 14);
    tree.insert(16, 16);
    tree.insert(15, 15);
    tree.insert(17, 17);
    tree.run(); 
    std::cout << std::endl;

    return (0);
}
