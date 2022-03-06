#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <utility>
#include <chrono>
#include <thread>
#include <cmath>

using namespace std;



class Timer {
private:
    chrono::time_point<chrono::system_clock> starting_point{};
    chrono::time_point<chrono::system_clock> ending_point{};
    long double counter = 0.0;
    int counter_increment = 0;
public:
    Timer() = default;

    void start() {this->starting_point = chrono::system_clock::now();}

    void end() {this->ending_point = chrono::system_clock::now();}

    long double get_result() {
        chrono::duration<double> elapsed = ending_point - starting_point;
        return elapsed.count();
    }

    void operator ++ (int nothing) {
        this->counter += this->get_result();
        this->counter_increment++;
    }

    long double get_average() const {
        return this->counter / this->counter_increment;
    }
};


template <typename T>
vector<vector<T>> split_vector(vector<T> vec, int parts) {
    vector <vector<T> > bunches;
    int bunch_size = vec.size() / parts;
    for(size_t i = 0; i < vec.size(); i += bunch_size) {
        auto last = std::min(vec.size(), i + bunch_size);
        bunches.emplace_back(vec.begin() + i, vec.begin() + last);
    }
    return bunches;
}


template <typename T>
ostream & operator << (ostream & out, vector<T> x) {
    for (long long i = 0; i < x.size(); i++) {
        out << x[i] << "\t";
    }
    return out;
}


template <typename T>
class Matrix {
protected:
    vector<vector<T>> data;
    int rows{};
    int columns{};
public:
    Matrix() = default;

    Matrix(int _rows, int _cols, T k=0) {
        for (int i = 0; i < _rows; i++) {
            vector<T> temp(_cols);
            for (int j = 0; j < _cols; j++) {
                temp[j] = k;
            }
            data.push_back(temp);
        }
        rows = _rows;
        columns = _cols;
    }

    Matrix(vector<vector<T>> inp) {
        data = inp;
    }

    Matrix<T> copy() {
        Matrix<T> res(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                res.data[i][j] = data[i][j];
            }
        }
        return res;
    }

    void out() {
        for (int i = 0; i < data.size(); i++) {
            for (int j = 0; j < data[0].size(); j++) {
                cout << data[i][j] << " ";
            }
            cout << endl;
        }
    }

    void in() {
        for (int i = 0; i < data.size(); i++) {
            for (int j = 0; j < data[0].size(); j++) {
                cin >> data[i][j];
            }
        }
    }

    vector<T> & operator [](int row) {
        return data[row];
    }
};


template <typename T>
vector<T> solve_equation_bottom(Matrix<T> A, vector<T> b) {
    //Решение СЛАУ. А — нижнетреугольная матрица.
    vector<T> res(b.size());
    res[0] = b[0] / A[0][0];
    for (int i = 1; i < res.size(); i++) {
        T temp = 0;
        for (int j = i - 1; j >= 0; j--) {
            temp += res[j] * A[i][j];
        }
        res[i] = (b[i] - temp) / A[i][i];
    }
    return res;
}


template <typename T>
vector<T> solve_equation_top(Matrix<T> A, vector<T> b) {
    //Решение СЛАУ. А — верхнетреугольная матрица.
    vector<T> res(b.size());
    int size = res.size();
    res[size - 1] = b[size - 1] / A[size - 1][size - 1];
    for (int i = size - 2; i >= 0; i--) {
        T temp = 0;
        for (int j = i + 1; j <= size - 1; j++) {
            temp += res[j] * A[i][j];
        }
        res[i] = (b[i] - temp) / A[i][i];
    }
    return res;
}




template <typename T>
void LU(Matrix<T> A, Matrix<T> & L, Matrix<T> & U, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                L[i][j] = 1;
            }
        }
    }

    for(int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i <= j) {
                T S = 0;
                for (int k = 0; k < i; k++) {
                    S += L[i][k] * U[k][j];
                }
                U[i][j] = A[i][j] - S;
            } else {
                T S = 0;
                for (int k = 0; k < j; k++) {
                    S += L[i][k] * U[k][j];
                }
                L[i][j] = (A[i][j] - S) / U[j][j];
            }
        }
    }

}


class Point {
protected:
    double x = 0;
    double y = 0;
    string name;

public:
    explicit Point(double x = 0, double y = 0, string name = "Unknown") {
        this->x = x;
        this->y = y;
        this->name = std::move(name);
    }

    double get_distance(const Point& another) {
        return sqrt(pow(this->x - another.x, 2) + pow(this->y - another.y, 2));
    }

    double get_x() {
        return x;
    }

    virtual void info() {
        cout << "Point named " << name << " located at: " << x << " " << y << endl;
    }
};



class InfChargedLine : public Point {
private:
    double x_;
    double y_;
    double U;
    long double tau;

public:
    InfChargedLine(double x, double y, double x_, double y_, string name, double U, long double tau = 0) : Point(x, y, std::move(name)) {
        this->x_ = x_;
        this->y_ = y_;
        this->U = U;
        this->tau = tau;
    }

    InfChargedLine mirror() {
        InfChargedLine res(x, -y, x_, -y_, name + " mirrored", U, -tau);
        return res;
    }

    void set_charge(long double _tau) {
        this->tau = _tau;
    }

    long double get_charge() {
        return this->tau;
    }

    long double kfi(InfChargedLine another) {
        double pi = 3.1415926;
        long double eps0 = 8.85 * pow(10, -12);
        long double interm = sqrt(pow(this->x - another.x_, 2) + pow(this->y + another.y_, 2)) /
                        sqrt(pow(this->x - another.x_, 2) + pow(this->y - another.y_, 2));
        long double res = (1.0 / (2.0 * pi * eps0)) * log(interm);
        return res;
    }

    double get_U() {
        return U;
    }

    void info() {
        cout << "Line named " << name << " located at: " << x << " " << y << " charged with " << tau << endl;
    }
};


class InfChargedLinesSystem {

public:
    InfChargedLinesSystem() = default;

    explicit InfChargedLinesSystem(vector<InfChargedLine> lines) {
        this->lines = std::move(lines);
    }
    vector<InfChargedLinesSystem> split(int N) {
        vector<vector<InfChargedLine>> splitted = split_vector(lines, N);
        vector<InfChargedLinesSystem> res;
        for (int i = 0; i < splitted.size(); i++) {
            res[i] = InfChargedLinesSystem(splitted[i]);
        }
        return res;
    }

    void mirror() {
        int size = lines.size();
        for (int i = 0; i < size; i++) {
            lines.push_back(lines[i].mirror());
        }
    }

    Matrix<long double> form_matrix() {
        vector<vector<long double>> res;
        for (int i = 0; i < lines.size(); i++) {
            vector<long double> row;
            for (int j = 0; j < lines.size(); j++) {
                row.push_back(lines[i].kfi(lines[j]));
            }
            res.push_back(row);
        }
        Matrix<long double> a(res);
        return a;
    }

    vector<long double> form_U_vector() {
        vector<long double> Uv;
        for (int i = 0; i < lines.size(); i++) {
            Uv.push_back(lines[i].get_U());
        }
        return Uv;
    }

    void calc_charges() {
        //TODO: change vectors to matrix
        Matrix<long double> A = form_matrix();
        vector<long double> U_vec = form_U_vector();
        Matrix<long double> L(lines.size(), lines.size(), 0);
        Matrix<long double> U(lines.size(), lines.size(), 0);
        int n = lines.size();
        LU(A, L, U, lines.size());
        vector<long double> y = solve_equation_bottom(L, U_vec);
        vector<long double> TAU = solve_equation_top(U, y);
        for (int i = 0; i < lines.size(); i++) {
            lines[i].set_charge(TAU[i]);
        }
    }

    void info() {
        for (auto & line : lines) {
            line.info();
        }
    }

    vector<long double> calcE(vector<Point> points) {
        //TODO: исправить
        long double pi = 3.1415926;
        long double eps0 = 8.85 * pow(10, -12);
        vector<long double> E;
        for (auto & target : points) {
            long double Et = 0;
            for (auto & line : this->lines) {
                InfChargedLine mirrored_line = line.mirror();
                long double Eplus = line.get_charge() / (2 * pi * eps0 * line.get_distance(target));
                long double Eminus = line.get_charge() / (2 * pi * eps0 * mirrored_line.get_distance(target));
                Et += Eplus + Eminus;
            }
            E.push_back(Et);
        }
        return E;
    }

    vector<InfChargedLine> lines;
};


class Cabel {
protected:
    double x, y, R0, H0, U;
    string name;

public:
    Cabel(double R0, double H0, double U, string name="None", double x = 0) {
        this->x = x;
        this->y = H0;
        this->H0 = H0;
        this->U = U;
        this->R0 = R0;
        this->name = std::move(name);
    }

    InfChargedLinesSystem eqCableSplit(int N, double H) {
        H = this->R0 - H;
        vector<InfChargedLine> charges;
        double deltA = 2.0 * 3.1415926 / (double)N;
        double L = deltA * (this->R0 - H);
        for (int i = 0; i < N; i++) {
            double _x = (this->R0 - H) * cos(deltA * i);
            double _y = this->H0 + (this->R0 - H) * sin(deltA * i);
            double _xkt = this->R0 * cos(deltA * i);
            double _ykt = this->H0 + this->R0 * sin(deltA * i);
            string line_name = this->name + "_" + to_string(i);
            charges.emplace_back(_x, _y, _xkt, _ykt, line_name, this->U);
        }
        return InfChargedLinesSystem(charges);
    }
};


void foo_slau(Cabel cabel, int start, int end, int N, long double * tt) {
    InfChargedLinesSystem sys = cabel.eqCableSplit(N, 0.01);
    Timer slau;
    slau.start();
    sys.calc_charges();
    slau.end();
    *tt += slau.get_result();
    vector<Point> targets;
    for (int x = start; x <= end; x += 1) {
        double true_x = (double) x / 100.0;
        targets.emplace_back(true_x, 0.5);
    }
    vector<long double> E = sys.calcE(targets);
}


void foo(Cabel cabel, int start, int end, int N) {
    InfChargedLinesSystem sys = cabel.eqCableSplit(N, 0.01);
    sys.calc_charges();
    vector<Point> targets;
    for (int x = start; x <= end; x += 1) {
        double true_x = (double) x / 100.0;
        targets.emplace_back(true_x, 0.5);
    }
    vector<long double> E = sys.calcE(targets);
}


void foo_split_n(Cabel cabel, int start, int end, InfChargedLinesSystem sys) {
    sys.calc_charges();
    vector<Point> targets;
    for (int x = start; x <= end; x += 1) {
        double true_x = (double) x / 100.0;
        targets.emplace_back(true_x, 0.5);
    }
    vector<long double> E = sys.calcE(targets);
}


void check_thread_over_reg_slau_separated() {  //TODO: run with 2, 4, 8 threads
    for (int N = 50; N <= 500; N+=50) {
        for (int start = 20000; start <= 40000; start += 5000) {
            Timer regular_timer;
            long double thread_timer = 0.0;
            long double * thread_timer_pointer = &thread_timer;

            for (int j = 1; j <= 40; j++) {

                Cabel cabel(0.08, 0.8, 1000, "test_cabel");
                InfChargedLinesSystem sys = cabel.eqCableSplit(N, 0.01);
                regular_timer.start();
                sys.calc_charges();
                regular_timer.end();
                regular_timer++;

                vector<Point> targets;
                for (int x = -start; x <= start; x += 1) {
                    double true_x = (double) x / 100.0;
                    targets.emplace_back(true_x, 0.5);
                }
                vector<long double> E = sys.calcE(targets);

                /*for (int k = 0; k < 2; k++) {
                    thread f(foo, cabel, start);
                    f.join();
                }*/
                thread s1(foo_slau, cabel, 0, start, N, thread_timer_pointer);
                thread s2(foo_slau, cabel, 0, start, N, thread_timer_pointer);
                /*thread s3(foo, cabel, 0, start / 4, N);
                thread s4(foo, cabel, 0, start / 4, N);
                thread s5(foo, cabel, 0, start / 4, N);
                thread s6(foo, cabel, 0, start / 4, N);
                thread s7(foo, cabel, 0, start / 4, N);
                thread s8(foo, cabel, 0, start / 4, N);*/
                s1.join();
                s2.join();
                /*s3.join();
                s4.join();
                s5.join();
                s6.join();
                s7.join();
                s8.join();*/

            }

            cout << regular_timer.get_average() << "\t";
            cout << thread_timer / 40.0 << "\t";
            cout << start << "\t";
            cout << N << endl;
        }
    }
}


void check_thread_over_reg_attempts_split_points() {  //TODO: run with 2, 4, 8 threads
    for (int N = 50; N <= 500; N+=50) {
        for (int start = 20000; start <= 30000; start += 15000) {
            Timer regular_timer;
            Timer thread_timer;

            for (int j = 1; j <= 40; j++) {
                regular_timer.start();

                Cabel cabel(0.08, 0.8, 1000, "test_cabel");
                InfChargedLinesSystem sys = cabel.eqCableSplit(N, 0.01);
                sys.calc_charges();

                vector<Point> targets;
                for (int x = -start; x <= start; x += 1) {
                    double true_x = (double) x / 100.0;
                    targets.emplace_back(true_x, 0.5);
                }
                vector<long double> E = sys.calcE(targets);

                regular_timer.end();
                regular_timer++;

                thread_timer.start();

                /*for (int k = 0; k < 2; k++) {
                    thread f(foo, cabel, start);
                    f.join();
                }
                thread s1(foo, cabel, 0, start / 4, N);
                thread s2(foo, cabel, 0, start / 4, N);
                thread s3(foo, cabel, 0, start / 4, N);
                thread s4(foo, cabel, 0, start / 4, N);
                thread s5(foo, cabel, 0, start / 4, N);
                thread s6(foo, cabel, 0, start / 4, N);
                thread s7(foo, cabel, 0, start / 4, N);
                thread s8(foo, cabel, 0, start / 4, N);
                s1.join();
                s2.join();
                s3.join();
                s4.join();
                s5.join();
                s6.join();
                s7.join();
                s8.join();*/

                thread_timer.end();
                thread_timer++;

            }

            cout << regular_timer.get_average() << "\t";
            cout << thread_timer.get_average() << "\t";
            cout << start << "\t";
            cout << N << endl;
        }
    }
}


void check_thread_over_reg_attempts_split_N() {  //TODO: run with 2, 4, 8 threads
    for (int N = 50; N <= 500; N+=50) {
        for (int start = 20000; start <= 30000; start += 15000) {
            Timer regular_timer;
            Timer thread_timer;

            for (int j = 1; j <= 40; j++) {

                Cabel cabel_(0.08, 0.8, 1000, "test_cabel");
                InfChargedLinesSystem sys_ = cabel_.eqCableSplit(N, 0.01);
                vector<InfChargedLine> lines1;
                vector<InfChargedLine> lines2;

                for (int i = 0; i < sys_.lines.size(); i++) {
                    if (i < sys_.lines.size() / 2) {
                        lines1.push_back(sys_.lines[i]);
                    } else {
                        lines2.push_back(sys_.lines[i]);
                    }
                }

                InfChargedLinesSystem part1(lines1);
                InfChargedLinesSystem part2(lines2);

                regular_timer.start();

                Cabel cabel(0.08, 0.8, 1000, "test_cabel");
                InfChargedLinesSystem sys = cabel.eqCableSplit(N, 0.01);
                sys.calc_charges();

                vector<Point> targets;
                for (int x = -start; x <= start; x += 1) {
                    double true_x = (double) x / 100.0;
                    targets.emplace_back(true_x, 0.5);
                }
                vector<long double> E = sys.calcE(targets);

                regular_timer.end();
                regular_timer++;

                thread_timer.start();

                thread s1(foo_split_n, cabel, -start, start, part1);
                thread s2(foo_split_n, cabel, -start, start, part2);

                thread_timer.end();
                thread_timer++;

            }

            cout << regular_timer.get_average() << "\t";
            cout << thread_timer.get_average() << "\t";
            cout << start << "\t";
            cout << N << endl;
        }
    }
}


void solve(vector<Point> targets, InfChargedLinesSystem sys) {
    vector<long double> E = sys.calcE(std::move(targets));
    cout << E << endl;
}


void check_thread_reg_equality() {

    Cabel cabel(0.08, 0.8, 1000, "test_cabel");
    InfChargedLinesSystem sys = cabel.eqCableSplit(3, 0.01);
    sys.calc_charges();

    vector<Point> targets;
    for (int x = -20; x <= 20; x += 1) {
        double true_x = (double) x / 100.0;
        targets.emplace_back(true_x, 0.5);
    }

    vector< vector<Point> > targets_bunches = split_vector(targets, 2);
    //TODO: fix (maybe by hands)


    cout << "regular" << endl;
    solve(targets, sys);
    thread f(solve, targets_bunches[0], sys);
    thread s(solve, targets_bunches[1], sys);
    cout << "first thread" << endl;
    f.join();
    cout << "second thread" << endl;
    s.join();
}


int main() {
    //check_thread_reg_equality();
    //check_thread_over_reg_slau_separated();
    //check_thread_over_reg_attempts_split_points();
    check_thread_over_reg_attempts_split_N();
}
