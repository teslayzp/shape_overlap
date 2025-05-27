class ShapeManager {
    constructor() {
        this.canvas = document.getElementById('shapeCanvas');
        if (!this.canvas) {
            console.error('Could not find canvas element');
            return;
        }
        
        this.ctx = this.canvas.getContext('2d');
        if (!this.ctx) {
            console.error('Could not get canvas context');
            return;
        }

        this.shapes = [];
        this.selectedShape = null;
        this.isDragging = false;
        this.lastX = 0;
        this.lastY = 0;

        // Set canvas size
        this.resizeCanvas();
        
        // Debug log
        console.log('Canvas initialized:', {
            width: this.canvas.width,
            height: this.canvas.height
        });

        // Event listeners
        this.setupEventListeners();
    }

    setupEventListeners() {
        // Canvas events
        this.canvas.addEventListener('mousedown', (e) => this.handleMouseDown(e));
        this.canvas.addEventListener('mousemove', (e) => this.handleMouseMove(e));
        this.canvas.addEventListener('mouseup', () => this.handleMouseUp());

        // Shape controls
        const shapeTypeSelect = document.getElementById('shapeType');
        if (shapeTypeSelect) {
            shapeTypeSelect.addEventListener('change', () => this.toggleControls());
        }

        // Auto shape type controls
        const autoShapeType = document.getElementById('autoShapeType');
        if (autoShapeType) {
            autoShapeType.addEventListener('change', () => {
                const controls = document.getElementById('autoPolygonControls');
                if (controls) {
                    controls.style.display = autoShapeType.value === 'circle' ? 'none' : 'grid';
                }
            });
        }

        // Buttons
        const addButton = document.getElementById('addShape');
        if (addButton) {
            addButton.addEventListener('click', () => this.addShape());
        }

        const checkButton = document.getElementById('checkOverlap');
        if (checkButton) {
            checkButton.addEventListener('click', () => this.checkOverlaps());
        }

        const clearButton = document.getElementById('clearShapes');
        if (clearButton) {
            clearButton.addEventListener('click', () => this.clearShapes());
        }

        const generateButton = document.getElementById('generateShapes');
        if (generateButton) {
            generateButton.addEventListener('click', () => this.generateAndCheckShapes());
        }

        // Initial toggle of controls
        this.toggleControls();
    }

    resizeCanvas() {
        if (!this.canvas) return;
        
        const container = this.canvas.parentElement;
        if (!container) {
            console.error('Could not find canvas container');
            return;
        }

        this.canvas.width = container.clientWidth - 40;
        this.canvas.height = 600;
        
        console.log('Canvas resized:', {
            width: this.canvas.width,
            height: this.canvas.height,
            containerWidth: container.clientWidth
        });

        this.draw();
    }

    generateAndCheckShapes() {
        // Validate inputs first
        const errors = this.validateGeneratorInputs();
        if (errors.length > 0) {
            const results = document.getElementById('overlapResults');
            if (results) {
                results.innerHTML = `<p style="color: var(--error-color)">Errors:<br>${errors.join('<br>')}</p>`;
            }
            return;
        }

        // Clear existing shapes
        this.clearShapes();

        try {
            const numShapes = parseInt(document.getElementById('numShapes').value);
            const shapeType = document.getElementById('autoShapeType').value;

            // Generate shapes
            for (let i = 0; i < numShapes; i++) {
                const shape = this.generateRandomShape(shapeType);
                if (!shape) {
                    throw new Error('Failed to generate shape');
                }
                this.shapes.push(shape);
            }

            // Update UI
            this.updateShapesList();
            this.draw();

            // Check overlaps
            this.checkOverlaps();
        } catch (error) {
            console.error('Error generating shapes:', error);
            const results = document.getElementById('overlapResults');
            if (results) {
                results.innerHTML = `<p style="color: var(--error-color)">Error generating shapes: ${error.message}</p>`;
            }
        }
    }

    validateGeneratorInputs() {
        const errors = [];
        
        // Validate number of shapes
        const numShapes = parseInt(document.getElementById('numShapes').value);
        if (isNaN(numShapes) || numShapes < 2 || numShapes > 100) {
            errors.push('Number of shapes must be between 2 and 100');
        }

        const shapeType = document.getElementById('autoShapeType').value;
        if (shapeType !== 'circle') {
            // Validate sides range for both regular and irregular polygons
            const minSides = parseInt(document.getElementById('minSides').value);
            const maxSides = parseInt(document.getElementById('maxSides').value);
            
            if (isNaN(minSides) || minSides < 3 || minSides > 50) {
                errors.push('Minimum sides must be between 3 and 50');
            }
            if (isNaN(maxSides) || maxSides < 3 || maxSides > 50) {
                errors.push('Maximum sides must be between 3 and 50');
            }
            if (minSides > maxSides) {
                errors.push('Minimum sides cannot be greater than maximum sides');
            }
        }

        // Validate size range
        const minSize = parseInt(document.getElementById('minSize').value);
        const maxSize = parseInt(document.getElementById('maxSize').value);
        
        if (isNaN(minSize) || minSize < 10 || minSize > 100) {
            errors.push('Minimum size must be between 10 and 100');
        }
        if (isNaN(maxSize) || maxSize < 10 || maxSize > 100) {
            errors.push('Maximum size must be between 10 and 100');
        }
        if (minSize > maxSize) {
            errors.push('Minimum size cannot be greater than maximum size');
        }

        return errors;
    }

    generateRandomShape(type) {
        try {
            // Calculate safe boundaries for shape placement
            const padding = 50;
            const maxX = this.canvas.width - padding;
            const maxY = this.canvas.height - padding;
            const minX = padding;
            const minY = padding;

            // Generate random position within safe boundaries
            const x = Math.random() * (maxX - minX) + minX;
            const y = Math.random() * (maxY - minY) + minY;

            // Get size range from inputs
            const minSize = parseInt(document.getElementById('minSize').value);
            const maxSize = parseInt(document.getElementById('maxSize').value);
            
            if (isNaN(minSize) || isNaN(maxSize)) {
                throw new Error('Invalid size values');
            }
            
            const radius = minSize + Math.random() * (maxSize - minSize);

            if (type === 'mixed') {
                const types = ['circle', 'regular_polygon', 'irregular_polygon'];
                type = types[Math.floor(Math.random() * types.length)];
            }

            if (type === 'circle') {
                return {
                    type: 'circle',
                    x: x,
                    y: y,
                    radius: radius,
                    isOverlapping: false
                };
            } else if (type === 'regular_polygon') {
                // Get sides range from inputs
                const minSides = parseInt(document.getElementById('minSides').value);
                const maxSides = parseInt(document.getElementById('maxSides').value);
                
                if (isNaN(minSides) || isNaN(maxSides)) {
                    throw new Error('Invalid sides values');
                }
                
                const sides = minSides + Math.floor(Math.random() * (maxSides - minSides + 1));

                return {
                    type: 'regular_polygon',
                    x: x,
                    y: y,
                    sides: sides,
                    radius: radius,
                    rotation: Math.random() * 360,
                    isOverlapping: false
                };
            } else if (type === 'irregular_polygon') {
                // Generate random number of vertices
                const minSides = parseInt(document.getElementById('minSides').value);
                const maxSides = parseInt(document.getElementById('maxSides').value);
                const numVertices = minSides + Math.floor(Math.random() * (maxSides - minSides + 1));
                
                // Generate vertices with random variations
                const vertices = [];
                const angleStep = (2 * Math.PI) / numVertices;
                
                for (let i = 0; i < numVertices; i++) {
                    const angle = i * angleStep + (Math.random() * 0.5 - 0.25) * angleStep; // Add some angle variation
                    const distance = radius * (0.7 + Math.random() * 0.6); // Vary the distance from center
                    vertices.push([
                        distance * Math.cos(angle),
                        distance * Math.sin(angle)
                    ]);
                }

                return {
                    type: 'irregular_polygon',
                    x: x,
                    y: y,
                    vertices: vertices,
                    isOverlapping: false
                };
            }
        } catch (error) {
            console.error('Error in generateRandomShape:', error);
            return null;
        }
    }

    setupCanvas() {
        // Make canvas responsive
        const resizeCanvas = () => {
            const container = this.canvas.parentElement;
            this.canvas.width = container.clientWidth - 40;
            this.canvas.height = 600;
            this.draw();
        };
        window.addEventListener('resize', resizeCanvas);
        resizeCanvas();
    }

    getMousePos(e) {
        const rect = this.canvas.getBoundingClientRect();
        return {
            x: e.clientX - rect.left,
            y: e.clientY - rect.top
        };
    }

    handleMouseDown(e) {
        const pos = this.getMousePos(e);
        for (let i = this.shapes.length - 1; i >= 0; i--) {
            if (this.isPointInShape(pos, this.shapes[i])) {
                this.selectedShape = this.shapes[i];
                this.isDragging = true;
                this.lastX = pos.x;
                this.lastY = pos.y;
                break;
            }
        }
    }

    handleMouseMove(e) {
        if (!this.isDragging || !this.selectedShape) return;

        const pos = this.getMousePos(e);
        const dx = pos.x - this.lastX;
        const dy = pos.y - this.lastY;

        this.selectedShape.x += dx;
        this.selectedShape.y += dy;

        this.lastX = pos.x;
        this.lastY = pos.y;

        this.draw();
    }

    handleMouseUp() {
        this.isDragging = false;
        this.selectedShape = null;
    }

    isPointInShape(point, shape) {
        if (shape.type === 'circle') {
            const dx = point.x - shape.x;
            const dy = point.y - shape.y;
            return Math.sqrt(dx * dx + dy * dy) <= shape.radius;
        }
        // For polygons, implement point-in-polygon check
        return false; // Temporary
    }

    addShape() {
        const type = document.getElementById('shapeType').value;
        const shape = {
            type,
            x: this.canvas.width / 2,
            y: this.canvas.height / 2
        };

        if (type === 'circle') {
            shape.radius = parseInt(document.getElementById('circleRadius').value);
        } else if (type === 'regular_polygon') {
            shape.sides = parseInt(document.getElementById('polygonSides').value);
            shape.radius = parseInt(document.getElementById('polygonRadius').value);
            shape.rotation = parseInt(document.getElementById('polygonRotation').value);
        }

        this.shapes.push(shape);
        this.updateShapesList();
        this.draw();
    }

    clearShapes() {
        this.shapes = [];
        this.updateShapesList();
        this.draw();
        const results = document.getElementById('overlapResults');
        if (results) {
            results.innerHTML = '';
        }
    }

    updateShapesList() {
        const list = document.getElementById('shapesList');
        if (!list) return;

        list.innerHTML = '';
        this.shapes.forEach((shape, index) => {
            const li = document.createElement('li');
            li.textContent = `${shape.type} ${index + 1}`;
            if (shape.isOverlapping) {
                li.style.color = 'var(--success-color)';
            }
            const deleteBtn = document.createElement('button');
            deleteBtn.textContent = 'Delete';
            deleteBtn.onclick = () => {
                this.shapes.splice(index, 1);
                this.updateShapesList();
                this.draw();
            };
            li.appendChild(deleteBtn);
            list.appendChild(li);
        });
    }

    checkOverlaps() {
        try {
            const results = document.getElementById('overlapResults');
            if (!results) return;

            results.innerHTML = '';
            let overlappingPairs = [];

            // Reset previous overlap states
            this.shapes.forEach(shape => shape.isOverlapping = false);

            // Check each pair of shapes
            for (let i = 0; i < this.shapes.length; i++) {
                for (let j = i + 1; j < this.shapes.length; j++) {
                    if (this.doShapesOverlap(this.shapes[i], this.shapes[j])) {
                        overlappingPairs.push([i, j]);
                        this.shapes[i].isOverlapping = true;
                        this.shapes[j].isOverlapping = true;
                    }
                }
            }

            // Display results
            if (overlappingPairs.length === 0) {
                results.innerHTML = '<p style="color: var(--error-color)">No overlapping shapes found.</p>';
            } else {
                let html = `<p style="color: var(--success-color)">Found ${overlappingPairs.length} overlapping pair(s):</p><ul>`;
                overlappingPairs.forEach(([i, j]) => {
                    html += `<li>Shape ${i + 1} overlaps with Shape ${j + 1}</li>`;
                });
                html += '</ul>';
                results.innerHTML = html;
            }

            // Redraw with highlights
            this.draw();
        } catch (error) {
            console.error('Error in checkOverlaps:', error);
            const results = document.getElementById('overlapResults');
            if (results) {
                results.innerHTML = `<p style="color: var(--error-color)">Error checking overlaps: ${error.message}</p>`;
            }
        }
    }

    doShapesOverlap(shape1, shape2) {
        // Circle-Circle overlap
        if (shape1.type === 'circle' && shape2.type === 'circle') {
            const dx = shape2.x - shape1.x;
            const dy = shape2.y - shape1.y;
            const distance = Math.sqrt(dx * dx + dy * dy);
            return distance < (shape1.radius + shape2.radius);
        }

        // Regular Polygon-Regular Polygon overlap
        if (shape1.type === 'regular_polygon' && shape2.type === 'regular_polygon') {
            return this.checkPolygonOverlap(
                this.getPolygonPoints(shape1),
                this.getPolygonPoints(shape2)
            );
        }

        // Circle-Regular Polygon overlap
        if (shape1.type === 'circle' && shape2.type === 'regular_polygon') {
            return this.checkCirclePolygonOverlap(shape1, shape2);
        }
        if (shape1.type === 'regular_polygon' && shape2.type === 'circle') {
            return this.checkCirclePolygonOverlap(shape2, shape1);
        }

        return false;
    }

    getPolygonPoints(polygon) {
        const points = [];
        const angleStep = (Math.PI * 2) / polygon.sides;
        const rotationRad = (polygon.rotation * Math.PI) / 180;

        for (let i = 0; i < polygon.sides; i++) {
            const angle = i * angleStep + rotationRad;
            points.push({
                x: polygon.x + polygon.radius * Math.cos(angle),
                y: polygon.y + polygon.radius * Math.sin(angle)
            });
        }
        return points;
    }

    checkPolygonOverlap(points1, points2) {
        // Check if any point from one polygon is inside the other
        for (let point of points1) {
            if (this.isPointInPolygon(point, points2)) return true;
        }
        for (let point of points2) {
            if (this.isPointInPolygon(point, points1)) return true;
        }

        // Check if any lines intersect
        for (let i = 0; i < points1.length; i++) {
            const i2 = (i + 1) % points1.length;
            for (let j = 0; j < points2.length; j++) {
                const j2 = (j + 1) % points2.length;
                if (this.doLinesIntersect(
                    points1[i], points1[i2],
                    points2[j], points2[j2]
                )) return true;
            }
        }

        return false;
    }

    checkCirclePolygonOverlap(circle, polygon) {
        const points = this.getPolygonPoints(polygon);
        
        // Check if circle center is inside polygon
        if (this.isPointInPolygon({x: circle.x, y: circle.y}, points)) return true;

        // Check if any polygon point is inside circle
        for (let point of points) {
            const dx = point.x - circle.x;
            const dy = point.y - circle.y;
            if (Math.sqrt(dx * dx + dy * dy) <= circle.radius) return true;
        }

        // Check if circle intersects with any polygon edge
        for (let i = 0; i < points.length; i++) {
            const j = (i + 1) % points.length;
            if (this.distanceToLineSegment(circle, points[i], points[j]) <= circle.radius) return true;
        }

        return false;
    }

    isPointInPolygon(point, polygonPoints) {
        let inside = false;
        for (let i = 0, j = polygonPoints.length - 1; i < polygonPoints.length; j = i++) {
            const xi = polygonPoints[i].x, yi = polygonPoints[i].y;
            const xj = polygonPoints[j].x, yj = polygonPoints[j].y;
            
            if (((yi > point.y) !== (yj > point.y)) &&
                (point.x < (xj - xi) * (point.y - yi) / (yj - yi) + xi)) {
                inside = !inside;
            }
        }
        return inside;
    }

    doLinesIntersect(p1, p2, p3, p4) {
        const denominator = (p4.y - p3.y) * (p2.x - p1.x) - (p4.x - p3.x) * (p2.y - p1.y);
        if (denominator === 0) return false;

        const ua = ((p4.x - p3.x) * (p1.y - p3.y) - (p4.y - p3.y) * (p1.x - p3.x)) / denominator;
        const ub = ((p2.x - p1.x) * (p1.y - p3.y) - (p2.y - p1.y) * (p1.x - p3.x)) / denominator;

        return ua >= 0 && ua <= 1 && ub >= 0 && ub <= 1;
    }

    distanceToLineSegment(circle, p1, p2) {
        const A = circle.x - p1.x;
        const B = circle.y - p1.y;
        const C = p2.x - p1.x;
        const D = p2.y - p1.y;

        const dot = A * C + B * D;
        const len_sq = C * C + D * D;
        let param = -1;

        if (len_sq !== 0) param = dot / len_sq;

        let xx, yy;

        if (param < 0) {
            xx = p1.x;
            yy = p1.y;
        } else if (param > 1) {
            xx = p2.x;
            yy = p2.y;
        } else {
            xx = p1.x + param * C;
            yy = p1.y + param * D;
        }

        const dx = circle.x - xx;
        const dy = circle.y - yy;
        return Math.sqrt(dx * dx + dy * dy);
    }

    draw() {
        if (!this.ctx || !this.canvas) return;
        
        // Clear canvas
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        
        // Draw shapes
        this.shapes.forEach(shape => this.drawShape(shape));
    }

    drawShape(shape) {
        this.ctx.beginPath();
        this.ctx.strokeStyle = shape.isOverlapping ? '#27ae60' : '#2c3e50';
        this.ctx.fillStyle = shape.isOverlapping ? 'rgba(39, 174, 96, 0.3)' : 'rgba(52, 152, 219, 0.5)';
        this.ctx.lineWidth = shape.isOverlapping ? 3 : 1;

        if (shape.type === 'circle') {
            this.ctx.arc(shape.x, shape.y, shape.radius, 0, Math.PI * 2);
        } else if (shape.type === 'regular_polygon') {
            const angleStep = (Math.PI * 2) / shape.sides;
            const rotationRad = (shape.rotation * Math.PI) / 180;
            
            for (let i = 0; i < shape.sides; i++) {
                const angle = i * angleStep + rotationRad;
                const x = shape.x + shape.radius * Math.cos(angle);
                const y = shape.y + shape.radius * Math.sin(angle);
                
                if (i === 0) {
                    this.ctx.moveTo(x, y);
                } else {
                    this.ctx.lineTo(x, y);
                }
            }
            this.ctx.closePath();
        } else if (shape.type === 'irregular_polygon') {
            const vertices = shape.vertices;
            if (vertices && vertices.length > 0) {
                this.ctx.moveTo(shape.x + vertices[0][0], shape.y + vertices[0][1]);
                for (let i = 1; i < vertices.length; i++) {
                    this.ctx.lineTo(shape.x + vertices[i][0], shape.y + vertices[i][1]);
                }
                this.ctx.closePath();
            }
        }

        this.ctx.fill();
        this.ctx.stroke();
    }

    toggleControls() {
        const shapeType = document.getElementById('shapeType').value;
        const circleControls = document.getElementById('circleControls');
        const polygonControls = document.getElementById('regularPolygonControls');

        if (circleControls && polygonControls) {
            circleControls.style.display = shapeType === 'circle' ? 'block' : 'none';
            polygonControls.style.display = shapeType === 'regular_polygon' ? 'block' : 'none';
        }
    }
}

// Initialize the shape manager when the page loads
window.addEventListener('load', () => {
    new ShapeManager();
}); 