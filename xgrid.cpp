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

    int num = 10;
    vector<vector<bool>> grid = set_init_matrix(num);
    int row = grid.size();
    int col = grid[0].size();

    auto is_valid = [&](int a, int b) {
        return (0 <= a) && (a < row) && (0 <= b) && (b <= col);
    };

    queue<pair<int, int>> nn, xt;
    vector<pair<int, int>> directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            for (pair<int, int> dir : directions) {
                int coordy = i + dir.first;
                int coordx = j + dir.second;
                if (is_valid(coordy, coordx) && grid[i][j] &&
                    !grid[coordx][coordy]) {
                    nn.push({coordy, coordx});
                }
            }
        }
    }
    do {
        while (!xt.empty()) {
            pair<int, int> coord = xt.front();
            xt.pop();
            for (pair<int, int> dir : directions) {
                int coordy = coord.first + dir.first;
                int coordx = coord.second + dir.second;
                if (is_valid(coordy, coordx) &&
                    grid[coord.first][coord.second] && !grid[coordy][coordx]) {
                    nn.push({coordx, coordy});
                }
            }
        }
        while (!nn.empty()) {
            if ((double)rand() / RAND_MAX) {
                xt.push(nn.front());
                grid[nn.front().first][nn.front().second] = true;
            }
            nn.pop();
        }
    } while (!xt.empty());

    // for (int i = 0; i < num; i++) {
    //     for (int j = 0; j < num; j++) {
    //         cout << grid[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    return 0;
}
