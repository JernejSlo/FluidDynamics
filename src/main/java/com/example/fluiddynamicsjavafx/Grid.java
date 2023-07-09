package com.example.fluiddynamicsjavafx;

import com.example.fluiddynamicsjavafx.Particle;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ConcurrentHashMap;

public class Grid {
    private ConcurrentHashMap<String, List<Particle>> grid;
    private int size;
    private double range;
    private int nrange;

    public Grid(int size, double range) {
        this.grid = new ConcurrentHashMap<>();
        this.size = size;
        this.range = range;
        this.nrange = (int) Math.round(range / size);
    }

    public void Insert(Particle p){
        int[] id = getCellID(p);
        String key = Arrays.toString(id);
        if (this.grid.get(key) == null){
            this.grid.put(key,new ArrayList<Particle>());
        }
        //System.out.println(this.grid.get(key));
        try{
            this.grid.get(key).add(p);
        } catch (Exception e) {
            System.out.println(key+" "+p+", "+e);
        }
    }

    public ArrayList<Particle> GetNeighbors(Particle particle){
        ArrayList<Particle> neighbors = new ArrayList<>();
        int[] id = getCellID(particle);
        for (int i = id[0] - (nrange); i <= id[0] + (nrange); i++){
            for (int j = id[1] - (nrange); j <= id[1] + (nrange); j++){
                String key = Arrays.toString(new int[]{i, j});
                if(grid.containsKey(key)){
                    for(Particle p: grid.get(key)){
                        if(distance(p.coords_, particle.coords_) <= range){
                            if (p.id != particle.id) {
                                neighbors.add(p);
                            }
                        }
                    }
                }
            }
        }
        return neighbors;
    }

    // utility method to calculate the Euclidean distance between two points
    private double distance(double [] p1, double [] p2){
        return Math.sqrt(Math.pow(p1[0]-p2[0],2) + Math.pow(p1[1]-p2[1],2));
    }

    // utility method to get cell ID from a particle's coordinates
    private int[] getCellID(Particle particle) {
        double[] coords = particle.coords_;
        int id_x = (int) Math.floor(coords[0] / size);
        int id_y = (int) Math.floor(coords[1] / size);
        return new int[] {id_x, id_y};
    }

    public static void main(String[] args) {
        Grid g = new Grid(3,4.);
        Particle p1 = new Particle(new double[]{1,1},1);
        Particle p2 = new Particle(new double[]{4,3},2);
        Particle p3 = new Particle(new double[]{1,22},3);
        g.Insert(p1);
        g.Insert(p2);
        g.Insert(p3);

        System.out.println("Searching neighbors for particle with coordinates [1, 1] and id 1:");
        for (Particle p : g.GetNeighbors(p1)){
            System.out.println(p);
        }

        System.out.println("\nSearching neighbors for particle with coordinates [4, 4] and id 2:");
        for (Particle p : g.GetNeighbors(p2)){
            System.out.println(p);
        }

        System.out.println("\nSearching neighbors for particle with coordinates [1, 22] and id 3:");
        for (Particle p : g.GetNeighbors(p3)){
            System.out.println(p);
        }
    }

}
