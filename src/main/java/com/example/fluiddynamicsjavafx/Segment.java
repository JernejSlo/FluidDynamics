package com.example.fluiddynamicsjavafx;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BrokenBarrierException;

public class Segment extends Thread{
    List<Particle> particles;
    Grid grid;
    public Segment(List<Particle> particles, Grid grid){
        this.particles = particles;
        this.grid = grid;
    }

    public void run(){
        //System.out.println(this.getId() + "is checking a array of the size " + particles.size());
        for(Particle particle : particles){
            Fluid2D.neighbors[particle.id] = (Fluid2D.grid.GetNeighbors(particle));
        }


        for (int i = 0; i < particles.size(); i++) {
            Particle currParticle = particles.get(i);
            Fluid2D.calculateDensity(currParticle);
            if (currParticle.getDensity() == 0){
                currParticle.setDensity(0.0001);
            }
        }

        for (int i = 0; i < particles.size(); i++) {
            Particle currParticle = particles.get(i);
            Fluid2D.calculatePressure(currParticle);
        }



        for (int i = 0; i < particles.size(); i++) {
            Particle currParticle = particles.get(i);
            Fluid2D.calculateForces(currParticle);
        }
        for (int i = 0; i < particles.size(); i++) {
            Particle currParticle = particles.get(i);
            currParticle.setVelocity(Fluid2D.newForces[currParticle.id]);
        }


        List<Particle> surfaceParticles = Fluid2D.getSurfaceParticles(particles);

        List<double[]> normalVectors = Fluid2D.getNormalVectors(surfaceParticles);

        List<Double> curvatureValues = Fluid2D.findCurves(surfaceParticles, normalVectors);

        /*for (int i = 0;i < surfaceParticles.size();i++) {
            double[] force = Fluid2D.calculateSurfaceTensionForce( curvatureValues.get(i),normalVectors.get(i),0.02, surfaceParticles.get(i).getPressure());
            if (Fluid2D.add(surfaceParticles.get(i).getVelocity(),force)[1]+"" != "NaN") {
                surfaceParticles.get(i).setVelocity(Fluid2D.add(surfaceParticles.get(i).getVelocity(),force));
            }
        }

         */
        for (Particle currParticle : particles){
            currParticle.setVelocity(Fluid2D.newForces[currParticle.id]);
            currParticle.setCoords_(Fluid2D.add(currParticle.getCoords_(), Fluid2D.mul(currParticle.getVelocity(),Fluid2D.timeStep)));
        }
        for (int i = 0; i < particles.size(); i++) {
            checkIfBounced(particles.get(i));
        }


    }

    public static void checkIfBounced(Particle currparticle){
        double []coords = currparticle.getCoords_();
        double []velocity = currparticle.getVelocity();
        if (coords[0] <= Fluid2D.MIN_X){

            double [] newVelocity = {-velocity[0]*Fluid2D.COR,velocity[1]};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] {0, currparticle.getCoords_()[1]});
        }
        if (coords[0] >= Fluid2D.MAX_X){

            double [] newVelocity = {-velocity[0]*Fluid2D.COR,velocity[1]};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] {Fluid2D.MAX_X, currparticle.getCoords_()[1]});
        }
        if (coords[1] <= 0){
            //double [] newVelocity = {velocity[0],-velocity[1]};
            double [] newVelocity = {velocity[0],-velocity[1]*Fluid2D.COR};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] { currparticle.getCoords_()[0],0});
        }
        if (coords[1] >= Fluid2D.MAX_Y){

            double [] newVelocity = {velocity[0],-velocity[1]*Fluid2D.COR};
            //double [] newVelocity = {velocity[0],-velocity[1]};
            currparticle.setVelocity(newVelocity);

            currparticle.setCoords_(new double[] { currparticle.getCoords_()[0],Fluid2D.MAX_Y});
        }
    }
}
