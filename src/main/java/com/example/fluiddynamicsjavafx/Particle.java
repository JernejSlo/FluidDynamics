package com.example.fluiddynamicsjavafx;

import javafx.scene.paint.Color;

public class Particle {

    public double[] coords_;
    private Color color = Color.LIGHTSKYBLUE;
    public int id;
    private int radius = Fluid2D.PARTICLE_RADIUS;
    private double [] velocity;
    private double density;
    private double pressure;

    public Particle(double[] coords, int id ) {
        coords_ = coords;
        this.id = id;
    }

    public Color getColor() {
        return color;
    }

    public void setColor(Color color) {
        this.color = color;
    }

    double get(int index) {
        return coords_[index];
    }
    double distance(Particle p2) {
        double dist = 0;
        for (int i = 0; i < coords_.length; ++i) {
            double d = coords_[i] - p2.coords_[i];
            dist += d * d;
        }
        return dist;
    }
    public String toString() {
        StringBuilder s = new StringBuilder("(");
        for (int i = 0; i < coords_.length; ++i) {
            if (i > 0)
                s.append(", ");
            s.append(coords_[i]);
        }
        s.append(')');
        return s.toString() + " with id: " + id;
    }
    public void show(){
        System.out.println("For particle " + this.id +" the values are:");
        System.out.println("coordinates: ("+this.coords_[0]+", "+this.coords_[1]+")");
        System.out.println("velocity: ("+this.velocity[0]+", "+this.velocity[1]+")");
        System.out.println("density: "+this.density);
        System.out.println("pressure: "+this.pressure);
    }

    public double[] getVelocity() {
        return velocity;
    }

    public double getDensity() {
        return density;
    }

    public void setDensity(double density) {
        this.density = density;
    }

    public double [] getCoords_() {
        return coords_;
    }

    public void setCoords_(double [] coords_) {
        this.coords_ = coords_;
    }

    public double getPressure() {
        return pressure;
    }

    public void setPressure(double pressure) {
        this.pressure = pressure;
    }

    public void setVelocity(double[] velocity) {
        this.velocity = velocity;
    }

    public int getRadius() {
        return radius;
    }

    public void setRadius(int radius) {
        this.radius = radius;
    }
}
