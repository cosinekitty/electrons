/*
    electrons.ts  -  Don Cross  -  http://cosinekitty.com
*/

/// <reference path="jquery.d.ts" />

module Electrons {
    function RadiansFromDegrees(d:number):number {
        return d * (Math.PI / 180);
    }

    class Vector {
        private x:number;
        private y:number;
        private z:number;

        public static Zero:Vector = new Vector(0, 0, 0);

        constructor(_x:number, _y:number, _z:number) {
            this.x = _x;
            this.y = _y;
            this.z = _z;
        }

        public getX():number { return this.x; }
        public getY():number { return this.y; }
        public getZ():number { return this.z; }

        public static Distance(a:Vector, b:Vector):number {
            var dx:number = b.getX() - a.getX();
            var dy:number = b.getY() - a.getY();
            var dz:number = b.getZ() - a.getZ();
            return Math.sqrt(dx*dx + dy*dy + dz*dz);
        }

        public static Cross(a:Vector, b:Vector):Vector {
            return new Vector(
                a.y*b.z - a.z*b.y,
                a.z*b.x - a.x*b.z,
                a.x*b.y - a.y*b.x
            );
        }

        public static Dot(a:Vector, b:Vector):number {
            return a.x*b.x + a.y*b.y + a.z*b.z;
        }

        public dot(other:Vector):number {
            return Vector.Dot(this, other);
        }

        public UnitVector():Vector {
            var mag:number = this.abs();
            return new Vector(this.x/mag, this.y/mag, this.z/mag);
        }

        public absSquared():number {
            return Vector.Dot(this, this);
        }

        public abs():number {
            return Math.sqrt(this.absSquared());
        }

        public add(other:Vector):Vector {
            return new Vector(
                this.x + other.x,
                this.y + other.y,
                this.z + other.z
            );
        }

        public sub(other:Vector):Vector {
            return new Vector(
                this.x - other.x,
                this.y - other.y,
                this.z - other.z
            );
        }

        public neg():Vector {
            return new Vector(-this.x, -this.y, -this.z);
        }

        public mul(scalar:number) {
            return new Vector(scalar*this.x, scalar*this.y, scalar*this.z);
        }

        public RotateX(a:number, b:number):Vector {
            return new Vector(this.x, a*this.y - b*this.z, b*this.y + a*this.z);
        }

        public RotateY(a:number, b:number):Vector {
            return new Vector(a*this.x + b*this.z, this.y, a*this.z - b*this.x);
        }

        public RotateZ(a:number, b:number):Vector {
            return new Vector(a*this.x - b*this.y, b*this.x - a*this.y, this.z);
        }
    }

    class RotationMatrix {
        private constructor(private r:Vector, private s:Vector, private t:Vector) {
        }

        public static Unrotated:RotationMatrix = new RotationMatrix(
            new Vector(1, 0, 0),
            new Vector(0, 1, 0),
            new Vector(0, 0, 1)
        );

        public RotateX(radians:number):RotationMatrix {
            let a:number = Math.cos(radians);
            let b:number = Math.sin(radians);
            return new RotationMatrix(
                this.r.RotateX(a, b),
                this.s.RotateX(a, b),
                this.t.RotateX(a, b));
        }

        public RotateY(radians:number):RotationMatrix {
            let a:number = Math.cos(radians);
            let b:number = Math.sin(radians);
            return new RotationMatrix(
                this.r.RotateY(a, b),
                this.s.RotateY(a, b),
                this.t.RotateY(a, b));
        }

        public RotateZ(radians:number):RotationMatrix {
            let a:number = Math.cos(radians);
            let b:number = Math.sin(radians);
            return new RotationMatrix(
                this.r.RotateZ(a, b),
                this.s.RotateZ(a, b),
                this.t.RotateZ(a, b));
        }

        public Rotate(v:Vector):Vector {
            return new Vector(v.dot(this.r), v.dot(this.s), v.dot(this.t));
        }
    }

    class CameraCoords {
        private hor:number;
        private ver:number;

        constructor (h:number, v:number) {
            this.hor = h;
            this.ver = v;
        }

        public getHor():number { return this.hor; }
        public getVer():number { return this.ver; }
    }

    class Display {
        public constructor(
            private pixelsWide:number,          // number of pixels wide in canvas
            private pixelsHigh:number,          // number of pixels high in canvas
            private zoomFactor:number,          // large values zoom in, small values zoom out
            private parallaxDistance:number)    // for scaling effect of distance on perspective
        {
        }

        public Erase(context:CanvasRenderingContext2D):void {
            context.clearRect(0, 0, this.pixelsWide, this.pixelsHigh);
            //context.strokeRect(0, 0, this.pixelsWide, this.pixelsHigh);
        }

        private GetCameraCoords(point:Vector):CameraCoords {
            var scale:number = this.pixelsWide * this.zoomFactor / (this.parallaxDistance - point.getZ());
            var h:number = scale * point.getX();
            var v:number = scale * point.getY();
            return new CameraCoords(this.pixelsWide/2 + h, this.pixelsHigh/2 - v);
        }

        public DrawSphere(
            context:CanvasRenderingContext2D,
            center:Vector,
            radius:number,
            frontColor:string,
            backColor:string = null,
            zLimit:number = null):number
        {
            // NOTE: This isn't quite right. The actual projection of a sphere
            // onto the pinhole camera screen is an ellipse, not a circle.
            // This will matter only for very low zoom with very close spheres
            // that are far off center.

            let origin:CameraCoords = this.GetCameraCoords(center);
            let rho:number = radius / this.parallaxDistance;
            let xp:number = rho * Math.sqrt(this.parallaxDistance*this.parallaxDistance - radius*radius);
            let zp:number = rho * radius;
            let tangent:Vector = center.add(new Vector(xp, 0, zp));
            let edge:CameraCoords = this.GetCameraCoords(tangent);
            let cradius:number = Math.abs(edge.getHor() - origin.getHor());

            context.beginPath();
            context.arc(origin.getHor(), origin.getVer(), cradius, 0, 2*Math.PI, true);
            if ((backColor !== null) && (zLimit !== null) && (center.getZ() < zLimit)) {
                context.strokeStyle = backColor;
            } else {
                context.strokeStyle = frontColor;
            }
            context.lineWidth = 1;
            context.stroke();

            return tangent.getZ();   // the z-value beneath which an electron is "around the bend"
        }

        public DrawLine(
            context:CanvasRenderingContext2D,
            startpoint:Vector,
            endpoint:Vector,
            color:string):void
        {
            let startcam:CameraCoords = this.GetCameraCoords(startpoint);
            let endcam:CameraCoords = this.GetCameraCoords(endpoint);

            context.beginPath();
            context.moveTo(startcam.getHor(), startcam.getVer());
            context.lineTo(endcam.getHor(), endcam.getVer());
            context.strokeStyle = color;
            context.lineWidth = 1;
            context.stroke();
        }
    }

    class Particle {
        private force:Vector;
        private position:Vector;

        public constructor(position:Vector) {
            this.position = position.UnitVector();
        }

        public GetPosition():Vector {
            return this.position;
        }

        public ResetForce():void {
            this.force = Vector.Zero;
        }

        public AddForce(other:Particle):void {
            // Force of electrically charged particles:
            // F = k*q1*q2/r^2.
            // We simplify the problem as:  F = 1/r^2.
            // Force is along the direction of the line passing through both.
            let dp:Vector = this.position.sub(other.position);
            let forcemag:number = 1.0 / dp.absSquared();
            let force:Vector = dp.UnitVector().mul(forcemag);
            this.force = this.force.add(force);
            other.force = other.force.sub(force);
        }

        public Migrate(positionShift:Vector):void {
            // Update the position of the particle.
            // Move the particle in the direction of the force, but constrain
            // to the surface of the unit sphere.
            this.position = this.position.add(positionShift).UnitVector();
        }

        public TangentialForce():Vector {
            // Calculate the tangential component of the force along
            // the sphere's surface.
            // Force = Force_radial + Force_tangential
            // Force_tangential = Force - Force_radial
            // So calculate radial component using dot product and subtract
            // to get tangential component.
            let radialForce:Vector = this.position.mul(this.force.dot(this.position));
            return this.force.sub(radialForce);
        }

        public Rotate(rotmat:RotationMatrix):void {
            this.position = rotmat.Rotate(this.position);
        }
    }

    class Simulation {
        private particleList:Particle[];
        private sphereCenter:Vector = new Vector(0, 0, 0);
        private sphereRadius:number = 1.0;
        private xAxis:Vector = new Vector(1, 0, 0);
        private yAxis:Vector = new Vector(0, 1, 0);
        private zAxis:Vector = new Vector(0, 0, 1);

        public EnableConnectNearestNeighbors:boolean = true;

        public constructor() {
            this.particleList = [];
        }

        public InsertParticle(p:Particle):void {
            this.particleList.push(p);
        }

        public RemoveParticle():void {
            if (this.particleList.length > 1) {
                this.particleList.pop();
            }
        }

        public ParticleCount():number {
            return this.particleList.length;
        }

        public Update():void {
            for (let i=0; i < this.particleList.length; ++i) {
                this.particleList[i].ResetForce();
            }

            for (let i=0; i < this.particleList.length - 1; ++i) {
                for (let j=i+1; j < this.particleList.length; ++j) {
                    this.particleList[i].AddForce(this.particleList[j]);
                }
            }

            let tangentialForceList:Vector[] = [];
            let maxForceMag:number = null;
            for (let i:number=0; i < this.particleList.length; ++i) {
                let tf:Vector = this.particleList[i].TangentialForce();
                let forceMag:number = tf.abs();
                tangentialForceList.push(tf);
                if ((maxForceMag === null) || (maxForceMag < forceMag)) {
                    maxForceMag = forceMag;
                }
            }

            let dt:number = 0.08 / maxForceMag;
            if (dt > 0.005) {
                dt = 0.005;
            }

            for (let i:number=0; i < this.particleList.length; ++i) {
                this.particleList[i].Migrate(tangentialForceList[i].mul(dt));
            }
        }

        public Render(display:Display):void {
            let context:CanvasRenderingContext2D = canvas.getContext('2d');
            display.Erase(context);
            let zbend:number = display.DrawSphere(context, this.sphereCenter, this.sphereRadius, '#eee');
            for (let i:number = 0; i < this.particleList.length; ++i) {
                display.DrawSphere(context, this.particleList[i].GetPosition(), 0.01, '#000', '#aaa', zbend);
            }

            if (this.EnableConnectNearestNeighbors) {
                // Find the smallest distance between any two particles.
                let minDistance:number = null;
                for (let i:number = 0; i < this.particleList.length - 1; ++i) {
                    let ipos:Vector = this.particleList[i].GetPosition();
                    for (let j:number = i+1; j < this.particleList.length; ++j) {
                        let jpos:Vector = this.particleList[j].GetPosition();
                        let distance:number = Vector.Distance(ipos, jpos);
                        if ((minDistance === null) || (distance < minDistance)) {
                            minDistance = distance;
                        }
                    }
                }

                if (minDistance !== null) {
                    // Connect all pairs of particles whose distance is not much larger than the minimum.
                    let threshold:number = 1.02 * minDistance;
                    for (let i:number = 0; i < this.particleList.length - 1; ++i) {
                        let ipos:Vector = this.particleList[i].GetPosition();
                        for (let j:number = i+1; j < this.particleList.length; ++j) {
                            let jpos:Vector = this.particleList[j].GetPosition();
                            let distance:number = Vector.Distance(ipos, jpos);
                            if (distance <= threshold) {
                                let color:string = (ipos.getZ() > 0 && jpos.getZ() > 0) ? '#000' : '#aca';
                                display.DrawLine(context, ipos, jpos, color);
                            }
                        }
                    }
                }
            }
        }

        public Rotate(rotmat:RotationMatrix):void {
            for (let p of this.particleList) {
                p.Rotate(rotmat);
            }
        }
    }

    var ballCountDiv:JQuery;
    var canvas:HTMLCanvasElement;
    var sim:Simulation;
    var display:Display;
    const FrameDelayMillis:number = 30;
    const ZoomFactor:number = 7;
    const ParallaxDistance:number = 15.0;
    const MinParticleCount:number = 1;
    const MaxParticleCount:number = 200;
    const InitialParticleCount:number = 12;
    var spinner:RotationMatrix = RotationMatrix.Unrotated.RotateY(RadiansFromDegrees(0.15));
    var initialTilt:RotationMatrix = RotationMatrix.Unrotated.RotateX(RadiansFromDegrees(-15.0));

    function AnimationFrame():void {
        sim.Render(display);
        sim.Update();
        sim.Rotate(spinner);
        window.setTimeout(AnimationFrame, FrameDelayMillis);
    }

    function RandomUnitVector():Vector {
        let x:number = 2*Math.random() - 1;
        let y:number = 2*Math.random() - 1;
        let z:number = 2*Math.random() - 1;
        return new Vector(x, y, z).UnitVector();
    }

    $(document).ready(function(){
        ballCountDiv = $('#BallCountDiv');
        canvas = <HTMLCanvasElement> document.getElementById('SimCanvas');
        sim = new Simulation();
        for (let i:number = 0; i < InitialParticleCount; ++i) {
            sim.InsertParticle(new Particle(RandomUnitVector()));
        }
        ballCountDiv.text(sim.ParticleCount());
        display = new Display(canvas.width, canvas.height, ZoomFactor, ParallaxDistance);
        sim.Rotate(initialTilt);
        $('#IncrementButton').click(function(){
            if (sim.ParticleCount() < MaxParticleCount) {
                sim.InsertParticle(new Particle(RandomUnitVector()));
                ballCountDiv.text(sim.ParticleCount());
            }
        });
        $('#DecrementButton').click(function(){
            if (sim.ParticleCount() > MinParticleCount) {
                sim.RemoveParticle();
                ballCountDiv.text(sim.ParticleCount());
            }
        });
        AnimationFrame();
    });
}
