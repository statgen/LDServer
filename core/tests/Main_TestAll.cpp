#include <gmock/gmock.h>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);

    ofstream config_file("redis-connection.txt");
    for (int i = 1; i < argc; ++i) {
        config_file << argv[i] << std::endl;
    }
    config_file.close();

    return RUN_ALL_TESTS();
}