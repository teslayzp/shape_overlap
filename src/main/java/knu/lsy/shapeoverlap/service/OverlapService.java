package knu.lsy.shapeoverlap.service;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.springframework.stereotype.Service;

import knu.lsy.shapeoverlap.model.Shape;

@Service
public class OverlapService {

    /**
     * Find all overlapping pairs of shapes
     */
    public List<OverlapResult> findOverlaps(List<Shape> shapes) {
        List<OverlapResult> overlaps = new ArrayList<>();

        for (int i = 0; i < shapes.size(); i++) {
            for (int j = i + 1; j < shapes.size(); j++) {
                Shape shape1 = shapes.get(i);
                Shape shape2 = shapes.get(j);

                if (shape1.overlaps(shape2)) {
                    overlaps.add(new OverlapResult(i, j, shape1, shape2));
                }
            }
        }

        return overlaps;
    }

    /**
     * Group overlapping shapes using Union-Find algorithm
     * This solves the chain grouping problem where A overlaps B, B overlaps C,
     * so A, B, and C should be in the same group
     */
    public List<List<Integer>> groupOverlappingShapes(List<Shape> shapes) {
        int n = shapes.size();
        UnionFind uf = new UnionFind(n);

        // Find all overlaps and union them
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (shapes.get(i).overlaps(shapes.get(j))) {
                    uf.union(i, j);
                }
            }
        }

        // Group shapes by their root parent
        Map<Integer, List<Integer>> groups = new HashMap<>();
        for (int i = 0; i < n; i++) {
            int root = uf.find(i);
            groups.computeIfAbsent(root, k -> new ArrayList<>()).add(i);
        }

        return new ArrayList<>(groups.values());
    }

    /**
     * Count total number of overlapping pairs
     */
    public int countOverlaps(List<Shape> shapes) {
        return findOverlaps(shapes).size();
    }

    /**
     * Check if a specific shape overlaps with any other shape in the list
     */
    public boolean hasAnyOverlap(Shape targetShape, List<Shape> otherShapes) {
        return otherShapes.stream().anyMatch(targetShape::overlaps);
    }

    /**
     * Find all shapes that overlap with a specific shape
     */
    public List<Shape> findOverlappingShapes(Shape targetShape, List<Shape> allShapes) {
        return allShapes.stream()
                .filter(shape -> !shape.equals(targetShape) && targetShape.overlaps(shape))
                .collect(Collectors.toList());
    }

    /**
     * Calculate overlap statistics
     */
    public OverlapStatistics calculateStatistics(List<Shape> shapes) {
        List<OverlapResult> overlaps = findOverlaps(shapes);
        List<List<Integer>> groups = groupOverlappingShapes(shapes);

        int totalShapes = shapes.size();
        int overlappingPairs = overlaps.size();
        int overlappingGroups = groups.size();
        int shapesInGroups = groups.stream()
                .filter(group -> group.size() > 1)
                .mapToInt(List::size)
                .sum();
        int isolatedShapes = totalShapes - shapesInGroups;

        return new OverlapStatistics(
                totalShapes,
                overlappingPairs,
                overlappingGroups,
                shapesInGroups,
                isolatedShapes
        );
    }

    /**
     * Union-Find (Disjoint Set) data structure for grouping overlapping shapes
     */
    private static class UnionFind {
        private final int[] parent;
        private final int[] rank;

        public UnionFind(int n) {
            parent = new int[n];
            rank = new int[n];
            for (int i = 0; i < n; i++) {
                parent[i] = i;
                rank[i] = 0;
            }
        }

        public int find(int x) {
            if (parent[x] != x) {
                parent[x] = find(parent[x]); // Path compression
            }
            return parent[x];
        }

        public void union(int x, int y) {
            int rootX = find(x);
            int rootY = find(y);

            if (rootX != rootY) {
                // Union by rank
                if (rank[rootX] < rank[rootY]) {
                    parent[rootX] = rootY;
                } else if (rank[rootX] > rank[rootY]) {
                    parent[rootY] = rootX;
                } else {
                    parent[rootY] = rootX;
                    rank[rootX]++;
                }
            }
        }
    }

    /**
     * Result class for overlap detection
     */
    public static class OverlapResult {
        private final int index1;
        private final int index2;
        private final Shape shape1;
        private final Shape shape2;

        public OverlapResult(int index1, int index2, Shape shape1, Shape shape2) {
            this.index1 = index1;
            this.index2 = index2;
            this.shape1 = shape1;
            this.shape2 = shape2;
        }

        public int getIndex1() { return index1; }
        public int getIndex2() { return index2; }
        public Shape getShape1() { return shape1; }
        public Shape getShape2() { return shape2; }

        @Override
        public String toString() {
            return String.format("Overlap between shapes %d and %d: %s overlaps %s",
                    index1, index2, shape1.getShapeType(), shape2.getShapeType());
        }
    }

    /**
     * Statistics class for overlap analysis
     */
    public static class OverlapStatistics {
        private final int totalShapes;
        private final int overlappingPairs;
        private final int overlappingGroups;
        private final int shapesInGroups;
        private final int isolatedShapes;

        public OverlapStatistics(int totalShapes, int overlappingPairs, int overlappingGroups,
                                 int shapesInGroups, int isolatedShapes) {
            this.totalShapes = totalShapes;
            this.overlappingPairs = overlappingPairs;
            this.overlappingGroups = overlappingGroups;
            this.shapesInGroups = shapesInGroups;
            this.isolatedShapes = isolatedShapes;
        }

        public int getTotalShapes() { return totalShapes; }
        public int getOverlappingPairs() { return overlappingPairs; }
        public int getOverlappingGroups() { return overlappingGroups; }
        public int getShapesInGroups() { return shapesInGroups; }
        public int getIsolatedShapes() { return isolatedShapes; }

        @Override
        public String toString() {
            return String.format(
                    "OverlapStatistics{totalShapes=%d, overlappingPairs=%d, overlappingGroups=%d, shapesInGroups=%d, isolatedShapes=%d}",
                    totalShapes, overlappingPairs, overlappingGroups, shapesInGroups, isolatedShapes
            );
        }
    }
}












