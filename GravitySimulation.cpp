#include <math.h>
#include <fstream>
#include <iostream>

#include "Eigen/Core"

using namespace std;

const double G = 1; // 6.67 * pow(10, -11);
const double M = 1.988 * pow(10, 30);
const double m = 5.97 * pow(10, 24);
const double aph = 1.5210 * pow(10, 11);
const double sem = 1.4960 * pow(10, 11);

double solarScale = aph;
double globScale = 1;
double scale = solarScale;
double timestepB = 100;
double logscale = 1 / log10(scale / timestepB + 1);

double eps = 0.1 * scale;

const int numParticles = 2;

char DEFAULTFILE[] = "GravitySim.txt";

ofstream file;

typedef Eigen::Matrix<double, 3, 3> particle;
typedef Eigen::Vector3<double> vect;

class GravitySimulation
{
private:
    double time = 0;
    double maxTimeStep;
    double curTimeStep;
    Eigen::Matrix<double, 1, numParticles> particleMasses;
    Eigen::Matrix<particle, 1, numParticles> particles;
    int q = 0;
    double minDist = scale;
public:
    GravitySimulation(double mts, char * filename);
    ~GravitySimulation();
    void Interaction();
    void Update();
    void doTimestep(double tmin = 0);
    void AddParticle(double m, double x, double y, double z, 
                     double vx=0, double vy=0, double vz=0);
    double kinetEner();
    double graviEner();
    void Record();
};

GravitySimulation::GravitySimulation(double mts, char * filename)
{
    file.open(filename);
    this->maxTimeStep = mts;
    this->curTimeStep = mts;
}

GravitySimulation::~GravitySimulation(){
    file.close();
}

void GravitySimulation::Interaction() {
    particle A, B;
    for (int i=0; i < numParticles; i++) {
        for (int j = i+1; j < numParticles; j++) {
            A = this->particles[i];
            B = this->particles[j];
            vect d = A.row(0) - B.row(0);
            double dabs = sqrt(d.dot(d));
            if (dabs < this->minDist) {
                this->minDist = dabs;
            }
            vect F = (G / pow(pow(dabs, 2) + pow(eps, 2), 1.5)) * d;
            this->particles[i].row(2) -= F * this->particleMasses[j];
            this->particles[j].row(2) += F * this->particleMasses[i];
        }
    }
    // +0.001 prevents timestep->0 as mindist->0
    this->curTimeStep = this->maxTimeStep * (log10(this->minDist/timestepB + 1) * logscale + .001);
}

void GravitySimulation::Update() {
    particle updateMatrix;
    updateMatrix << 1, this->curTimeStep, 0, 
                    0, 1, this->curTimeStep,
                    0, 0, 0;
    this->particles = updateMatrix * this->particles;
    this->minDist = scale;
    this->time += this->curTimeStep;
}

void GravitySimulation::doTimestep(double tmin) {
    if (tmin == 0) {
        this->Interaction();
        this->Update();
    } else {
        double t = 0.0;
        while (t <= tmin) {
            this->Interaction();
            if ((t+this->curTimeStep) > tmin) {
                this->curTimeStep = tmin - t;
                this->Update();
                break;
            }
            this->Update();
            t += this->curTimeStep;
        }
    }
}

void GravitySimulation::AddParticle(double m, double x, double y, double z, 
                     double vx, double vy, double vz) {
    if (this->q >= numParticles) {
        throw invalid_argument("Too many particles added!");
    } else {
        this->particleMasses[this->q] = m;
        particle p;
        p << x, y, z, vx, vy, vz, 0, 0, 0;
        this->particles[this->q] = p;
        this->q++;
    }
}

double GravitySimulation::kinetEner() {
    double kt = 0;
    for (int i = 0; i < numParticles; i++) {
        kt += 0.5 * (particleMasses[i] * (particles[i].row(1).dot(particles[i].row(1))));
    }
    return kt;
}

double GravitySimulation::graviEner() {
    particle A, B;
    double Et=0;
    for (int i=0; i < numParticles; i++) {
        for (int j = i+1; j < numParticles; j++) {
            A = this->particles[i];
            B = this->particles[j];
            vect d = A.row(0) - B.row(0);
            // dont need to sqrt as squaring straight after?
            double dabs = sqrt(d.dot(d));
            Et -= (G * this->particleMasses[i] * this->particleMasses[j] / sqrt(pow(dabs, 2) + pow(eps, 2)));
        }
    }
    return Et;
}

void GravitySimulation::Record() {
    file << this->time << " " << this->kinetEner() << " " << this->graviEner() << " ";
    for (particle p : this->particles) {
        file << p.row(0) << " ";
        //cout << p.row(0) << " " << p.row(1) << " || ";
    }
    file << endl;
    //cout << endl;
}

void _2BodyExample() {
    GravitySimulation a = GravitySimulation(2, DEFAULTFILE);
    double vplan = sqrt(1.0 * 100.0 * (2.0/10.0 - 1.0/10.0));
    a.AddParticle(100.0, 0.0,   0.0, 0.0, 0.0, -m/M*vplan, 0.0);
    a.AddParticle(1.0, 10.0, 0.0, 0.0, 0.0, vplan, 0.0);
    for (int i = 0; i < 100; i++) {
        a.Record();
        a.doTimestep(1);
    }
    a.Record();
    //use python to animate and graph
}

main() {
    _2BodyExample();
}