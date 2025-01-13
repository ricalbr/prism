#include <bits/stdc++.h>

using namespace std;

// Calcolo combinatorio
long long comb(int n, int k) {
    if (k > n || k < 0)
        return 0;
    if (k == 0 || k == n)
        return 1;
    vector<long long> c(k + 1, 0);
    c[0] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = min(i, k); j > 0; --j) {
            c[j] = c[j] + c[j - 1];
        }
    }
    return c[k];
}

// Memoizzazione personalizzata per P
unordered_map<string, double> memo;

string make_key(int D, double eta, double dcr, int N, int C, int K) {
    ostringstream oss;
    oss << D << "," << eta << "," << dcr << "," << N << "," << C << "," << K;
    return oss.str();
}

double P(int D, double eta, double dcr, int N, int C, int K) {
    if (N < 0 || C < 0 || C > N)
        return 0.0;
    if (C == 0 && N == 0) {
        return comb(D, K) * pow(dcr, K) * pow(1 - dcr, D - K);
    }

    string key = make_key(D, eta, dcr, N, C, K);
    if (memo.find(key) != memo.end()) {
        return memo[key];
    }

    double term1 = ((1 - eta) + eta * (K + C) / D) * P(D, eta, dcr, N - 1, C, K);
    double term2 = (eta * (D - (K + C - 1)) / D) * P(D, eta, dcr, N - 1, C - 1, K);
    double result = term1 + term2;

    memo[key] = result;
    return result;
}

double analytic(int D, double eta, double dcr, int N, int L) {
    double sum = 0.0;
    for (int i = 0; i <= L; ++i) {
        sum += P(D, eta, dcr, N, i, L - i);
    }
    return sum;
}

vector<vector<double>> get_matrix(int num_det, double eta, double dcr, double xtk) {
    vector<vector<double>> V_th(num_det + 1, vector<double>(num_det + 1, 0.0));

    for (int n = 0; n <= num_det; ++n) {
        for (int m = 0; m <= num_det; ++m) {
            V_th[n][m] = analytic(num_det, eta, dcr, m, n);
        }
    }
    return V_th;
}

void export_matrix_to_file(const vector<vector<double>> &matrix, const string &filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Errore nell'apertura del file " << filename << endl;
        return;
    }

    file << setprecision(15) << scientific; // Set high precision with scientific notation
    for (const auto &row : matrix) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1)
                file << '\t'; // Aggiunge un tab se non è l'ultimo elemento
        }
        file << '\n';
    }

    file.close();
    cout << "Matrice esportata con successo nel file " << filename << endl;
}

int main() {
    int num_det = 64;
    double eta = 0.8;
    double dcr_min = 0.001, dcr_max = 0.01;
    double xtk = 0.0;

    std::vector<double> dcr(num_det);
    double log_max = std::log10(dcr_max), log_min = std::log10(dcr_min);
    double step = (log_max - log_min) / (num_det - 1);
    for (int i = 0; i < num_det; ++i) {
        dcr[i] = std::pow(10.0, log_min + i * step);
    }
    auto mean = std::accumulate(std::begin(dcr), std::end(dcr), 0.0) / std::size(dcr);
    auto V_th = get_matrix(num_det, eta, mean, xtk);

    export_matrix_to_file(V_th, "V.txt");

    return 0;
}
