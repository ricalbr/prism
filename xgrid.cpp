#include <bits/stdc++.h>

const double EVENT_PROB = 0.2;
const double XTALK_PROB = 0.1;

using namespace std;
typedef long long ll;

vector<vector<bool>> set_init_matrix(int num) {
    vector<vector<bool>> mat(num, vector<bool>(num));

    for (int i = 0; i < num; i++) {
        for (int j = 0; j < num; j++) {
            mat[i][j] = ((double)rand() / RAND_MAX) < EVENT_PROB;
        }
    }
    return mat;
}

int32_t main(int argc, char *argv[]) {

    int num = 10000;
    auto start_time = std::chrono::high_resolution_clock::now();
    vector<vector<bool>> grid = set_init_matrix(num);
    int row = grid.size();
    int col = grid[0].size();

    auto is_valid = [&](int a, int b) {
        return (0 <= a) && (a < row) && (0 <= b) && (b < col);
    };

    queue<pair<int, int>> nn, xt;
    vector<pair<int, int>> directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

    for (int y = 0; y < row; y++) {
        for (int x = 0; x < col; x++) {
            if (!grid[y][x])
                continue;
            for (auto &[dy, dx] : directions) {
                int ny = y + dy, nx = x + dx;
                if (is_valid(ny, nx) && !grid[ny][nx]) {
                    nn.push({ny, nx});
                }
            }
        }
    }
    do {
        while (!xt.empty()) {
            auto [y, x] = xt.front();
            xt.pop();
            if (grid[y][x])
                continue;
            for (auto &[dy, dx] : directions) {
                int ny = y + dy, nx = x + dx;
                if (is_valid(ny, nx) && !grid[ny][nx]) {
                    nn.push({ny, nx});
                }
            }
        }
        while (!nn.empty()) {
            if ((double)rand() / RAND_MAX < XTALK_PROB) {
                xt.push(nn.front());
                grid[nn.front().first][nn.front().second] = true;
            }
            nn.pop();
        }
    } while (!xt.empty());

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << " s\n";

    // for (int i = 0; i < num; i++) {
    //     for (int j = 0; j < num; j++) {
    //         cout << grid[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    return 0;
}
