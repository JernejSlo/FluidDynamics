package com.example.fluiddynamicsjavafx;

// The Particle2D class represents a single fluid particle
class Particle {
    private double [] position;
    private double [] velocity;
    private double density;
    private double pressure;

    public Particle(double x, double y, double vx, double vy) {
        position = new double [2];
        position[0] = x;
        position[1] = y;
        velocity = new double [2];
        velocity[0] = vx;
        velocity[1] = vy;
    }



    public double [] getPosition() {
        return position;
    }

    public double [] getVelocity() {
        return velocity;
    }

    public double getDensity() {
        return density;
    }

    public double getPressure() {
        return pressure;
    }

    public void setPosition(double [] pos) {
        position = pos;
    }

    public void setVelocity(double [] vel) {
        velocity = vel;
    }

    public void setDensity(double dens) {
        density = dens;
    }

    public void setPressure(double press) {
        pressure = press;
    }
}

