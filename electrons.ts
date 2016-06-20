/*
    electrons.ts  -  Don Cross  -  http://cosinekitty.com
*/

module Electrons {
    'use strict';

    function RadiansFromDegrees(d:number):number {
        return d * (Math.PI / 180);
    }

    class Vector {
        readonly x:number;
        readonly y:number;
        readonly z:number;

        public static Zero:Vector = new Vector(0, 0, 0);

        constructor(_x:number, _y:number, _z:number) {
            this.x = _x;
            this.y = _y;
            this.z = _z;
        }

        public static Distance(a:Vector, b:Vector):number {
            var dx:number = b.x - a.x;
            var dy:number = b.y - a.y;
            var dz:number = b.z - a.z;
            return Math.sqrt(dx*dx + dy*dy + dz*dz);
        }

        public static Midpoint(a:Vector, b:Vector):Vector {
            return new Vector((a.x + b.x)/2, (a.y + b.y)/2, (a.z + b.z)/2);
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
        readonly hor:number;
        readonly ver:number;

        constructor (h:number, v:number) {
            this.hor = h;
            this.ver = v;
        }
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
            var scale:number = this.pixelsWide * this.zoomFactor / (this.parallaxDistance - point.z);
            var h:number = scale * point.x;
            var v:number = scale * point.y;
            return new CameraCoords(this.pixelsWide/2 + h, this.pixelsHigh/2 - v);
        }

        public DrawSphere(
            context:CanvasRenderingContext2D,
            center:Vector,
            radius:number,
            color:string):number
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
            let cradius:number = Math.abs(edge.hor - origin.hor);

            context.beginPath();
            context.arc(origin.hor, origin.ver, cradius, 0, 2*Math.PI, true);
            context.strokeStyle = color;
            context.lineWidth = 1;
            context.stroke();

            return tangent.z;   // the z-value beneath which an electron is "around the bend"
        }

        public DrawLine(
            context:CanvasRenderingContext2D,
            startpoint:Vector,
            endpoint:Vector,
            startcolor:string,
            endcolor:string,
            linedash:number[]):void
        {
            let startcam:CameraCoords = this.GetCameraCoords(startpoint);
            let endcam:CameraCoords = this.GetCameraCoords(endpoint);

            let gradient:CanvasGradient = context.createLinearGradient(
                startcam.hor, startcam.ver,
                endcam.hor, endcam.ver);

            gradient.addColorStop(0, startcolor);
            gradient.addColorStop(1, endcolor);

            context.setLineDash(linedash);
            context.beginPath();
            context.moveTo(startcam.hor, startcam.ver);
            context.lineTo(endcam.hor, endcam.ver);
            context.strokeStyle = gradient;
            context.lineWidth = 1;
            context.stroke();
            context.setLineDash([]);
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
            let n:number = this.particleList.length;
            for (let i=0; i < n; ++i) {
                this.particleList[i].ResetForce();
            }

            for (let i=0; i < n-1; ++i) {
                for (let j=i+1; j < n; ++j) {
                    this.particleList[i].AddForce(this.particleList[j]);
                }
            }

            let tangentialForceList:Vector[] = [];
            let maxForceMag:number = 0;
            for (let i:number=0; i < n; ++i) {
                let tf:Vector = this.particleList[i].TangentialForce();
                let forceMag:number = tf.abs();
                tangentialForceList.push(tf);
                maxForceMag = Math.max(maxForceMag, forceMag);
            }

            if (maxForceMag > 0) {
                // We want to move each electron a small distance compared width
                // the average distance between particles.
                // As the number of particles on the sphere increases, the
                // average distance between them goes down as an inverse square root.
                let dt:number = 0.005 / (maxForceMag * Math.sqrt(n));
                for (let i:number=0; i < n; ++i) {
                    this.particleList[i].Migrate(tangentialForceList[i].mul(dt));
                }
            }
        }

/*
        private DrawArrow(
            context:CanvasRenderingContext2D,
            x1:number, y1:number, x2:number, y2:number):void
        {
            context.beginPath();
            context.moveTo(x1, y1);
            context.lineTo(x2, y1);
            context.lineTo((x1+x2)/2, y2);
            context.lineTo(x1, y1);
            context.strokeStyle = '#000';
            context.lineWidth = 1;
            context.stroke();
        }

        private CountControlX:number = 730;
        private CountControlY:number = 700;
        private UpArrowX1 = this.CountControlX + 10;
        private UpArrowX2 = this.CountControlX + 40;
        private UpArrowY1 = this.CountControlY - 80;
        private UpArrowY2 = this.CountControlY - 50;

        private DownArrowX1 = this.CountControlX + 10;
        private DownArrowX2 = this.CountControlX + 40;
        private DownArrowY1 = this.CountControlY + 40;
        private DownArrowY2 = this.CountControlY + 70;

        private DrawParticleCountControls(context:CanvasRenderingContext2D):void {
            context.font = '20px sans';
            context.fillText('n=' + this.ParticleCount(), this.CountControlX, this.CountControlY);

            this.DrawArrow(context, this.UpArrowX1, this.UpArrowY2, this.UpArrowX2, this.UpArrowY1);
            this.DrawArrow(context, this.DownArrowX2, this.DownArrowY1, this.DownArrowX1, this.DownArrowY2);
        }
*/

        public CanvasMouseClick(x:number, y:number):void {
/*
            if ((x >= this.UpArrowX1) && (x <= this.UpArrowX2) && (y >= this.UpArrowY1) && (y <= this.UpArrowY2)) {
                this.InsertParticle(new Particle(RandomUnitVector()));
            }

            if ((x >= this.DownArrowX1) && (x <= this.DownArrowX2) && (y >= this.DownArrowY1) && (y <= this.DownArrowY2)) {
                this.RemoveParticle();
            }
*/
        }

        public AdjustParticleCount(newCount:number):void {
            while (this.ParticleCount() < newCount) {
                this.InsertParticle(new Particle(RandomUnitVector()));
            }

            while (this.ParticleCount() > newCount) {
                this.RemoveParticle();
            }
        }

        public Render(display:Display):void {
            let context:CanvasRenderingContext2D = canvas.getContext('2d');
            display.Erase(context);
            //this.DrawParticleCountControls(context);

            let zbend:number = display.DrawSphere(context, this.sphereCenter, this.sphereRadius, '#eee');
            for (let p of this.particleList) {
                display.DrawSphere(context, p.GetPosition(), 0.01, this.PointColor(p.GetPosition(), zbend));
            }

            let isConnectedIndex:{[index:number]:boolean} = {};
            let linedash = [[], [1,3], [1,6], [1,9], [1,13]];
            for (let level=0; level < 5; ++level) {
                // Find the smallest distance between any two (as-yet unconnected) particles.
                let nextIsConnectedIndex:{[index:number]:boolean} = {};
                let minDistance:number = null;
                for (let i:number = 0; i < this.particleList.length - 1; ++i) {
                    if (!isConnectedIndex[i]) {
                        let ipos:Vector = this.particleList[i].GetPosition();
                        for (let j:number = i+1; j < this.particleList.length; ++j) {
                            if (!isConnectedIndex[j]) {
                                let jpos:Vector = this.particleList[j].GetPosition();
                                let distance:number = Vector.Distance(ipos, jpos);
                                if ((minDistance === null) || (distance < minDistance)) {
                                    minDistance = distance;
                                }
                            }
                        }
                    }
                }

                if (minDistance === null) {
                    break;      // nothing more to connect
                }

                // Connect all pairs of particles whose distance is not much larger than the minimum.
                let threshold:number = 1.01 * minDistance;
                for (let i:number = 0; i < this.particleList.length - 1; ++i) {
                    if (!isConnectedIndex[i]) {
                        let ipos:Vector = this.particleList[i].GetPosition();
                        let icolor:string = this.PointColor(ipos, zbend);
                        for (let j:number = i+1; j < this.particleList.length; ++j) {
                            if (!isConnectedIndex[j]) {
                                let jpos:Vector = this.particleList[j].GetPosition();
                                let distance:number = Vector.Distance(ipos, jpos);
                                if (distance <= threshold) {
                                    let jcolor:string = this.PointColor(jpos, zbend);
                                    display.DrawLine(context, ipos, jpos, icolor, jcolor, linedash[level]);
                                    nextIsConnectedIndex[i] = nextIsConnectedIndex[j] = true;
                                }
                            }
                        }
                    }
                }

                for (let i in nextIsConnectedIndex) {
                    isConnectedIndex[i] = true;
                }
            }
        }

        private ColorRound(x:number):number {
            return Math.round(255 * Math.min(1, Math.max(0, x)));
        }

        private PointColor(p:Vector, zlimit:number):string {
            let frac:number = (zlimit - p.z) / (zlimit - (-1.1));
            let red:number = this.ColorRound(0.2 + (0.8*frac));
            let blue:number = this.ColorRound(1.2*frac);
            let green:number = this.ColorRound(frac);
            return 'rgb(' + red + ',' + green + ',' + blue + ')';
        }

        public Rotate(rotmat:RotationMatrix):void {
            for (let p of this.particleList) {
                p.Rotate(rotmat);
            }
        }
    }

    var canvas:HTMLCanvasElement;
    var sim:Simulation;
    var display:Display;
    const UpdatesPerFrame:number = 10;
    const FrameDelayMillis:number = 30;
    const ZoomFactor:number = 7;
    const ParallaxDistance:number = 15.0;
    const MinParticleCount:number = 1;
    const MaxParticleCount:number = 200;
    const InitialParticleCount:number = 14;
    var ySpinner:RotationMatrix = RotationMatrix.Unrotated.RotateY(RadiansFromDegrees(0.15));
    var xSpinner:RotationMatrix = RotationMatrix.Unrotated.RotateX(RadiansFromDegrees(0.0377));
    var initialTilt:RotationMatrix = RotationMatrix.Unrotated.RotateX(RadiansFromDegrees(-15.0));

    function AnimationFrame():void {
        sim.Render(display);
        for (let i=0; i < UpdatesPerFrame; ++i) {
            sim.Update();
        }
        sim.Rotate(ySpinner);
        sim.Rotate(xSpinner);
        window.setTimeout(AnimationFrame, FrameDelayMillis);
    }

    function RandomUnitVector():Vector {
        let x:number = 2*Math.random() - 1;
        let y:number = 2*Math.random() - 1;
        let z:number = 2*Math.random() - 1;
        return new Vector(x, y, z).UnitVector();
    }

    function OnCanvasClick(ev:MouseEvent) {
        let x:number = ev.pageX - canvas.offsetLeft;
        let y:number = ev.pageY - canvas.offsetTop;
        sim.CanvasMouseClick(x, y);
    }

    function OnEditParticleCount() {
        var particleCountEdit = <HTMLInputElement> document.getElementById('ParticleCountEditBox');
        var errorMessageDiv = document.getElementById('ErrorMessageDiv');
        var text = particleCountEdit.value;
        if (text.match(/^[0-9]{1,4}$/)) {
            var count = parseInt(text);
            if (count >= MinParticleCount && count <= MaxParticleCount) {
                sim.AdjustParticleCount(count);
                errorMessageDiv.textContent = '';
                particleCountEdit.blur();
                return;
            }
        }
        errorMessageDiv.textContent = 'Invalid number of particles. Must be an integer in the range ' +
            MinParticleCount + ' to ' + MaxParticleCount;
    }

    window.onload = function() {
        canvas = <HTMLCanvasElement> document.getElementById('SimCanvas');
        canvas.addEventListener('click', OnCanvasClick, false);
        sim = new Simulation();
        for (let i:number = 0; i < InitialParticleCount; ++i) {
            sim.InsertParticle(new Particle(RandomUnitVector()));
        }

        var particleCountEdit = <HTMLInputElement> document.getElementById('ParticleCountEditBox');
        particleCountEdit.value = InitialParticleCount.toFixed();
        particleCountEdit.onblur = OnEditParticleCount;
        particleCountEdit.setAttribute('min', MinParticleCount.toFixed());
        particleCountEdit.setAttribute('max', MaxParticleCount.toFixed());

        particleCountEdit.onkeypress = function(evt) {
            if (evt.keyCode === 13) {
                OnEditParticleCount();
                event.preventDefault();
                return false;
            }
            return true;
        }

        display = new Display(canvas.width, canvas.height, ZoomFactor, ParallaxDistance);
        sim.Rotate(initialTilt);
        AnimationFrame();
    }
}
